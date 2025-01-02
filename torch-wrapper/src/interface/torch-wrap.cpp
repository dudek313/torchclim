#include <onnxruntime_cxx_api.h>

#include <iostream>
#include <memory>
#include <cassert>
#include <array>
#include <sstream>

// helper function prototyping so it doesn't clutter the file
extern "C" void print_tensor(Ort::Value& tensor);
extern "C" std::string get_type_str(ONNXTensorElementDataType dtype);

// not sure if this is needed/defined again in torch_wrap_predict
static thread_local Ort::Session session(nullptr);
// can leave the below alone, state storage
static thread_local bool is_model_loaded = false;

enum InitState {
	MODEL_LOADING_FAILED = -1,
	MODEL_LOADED_SUCCESSFULLY = 1,
	MODEL_ALREADY_LOADED = 2
};

// this function is complete, no real modifications other than the session definition perhaps
extern "C"
int torch_wrap_create(const char* script_path) {
	if( is_model_loaded ) {
		return MODEL_ALREADY_LOADED;
	}

	assert( script_path != NULL );

	try {
		// Deserialize the ScriptModule from a file using Ort::Session.
		Ort::Env env;
		session = Ort::Session(env, script_path, Ort::SessionOptions{ nullptr });
		is_model_loaded = true;
	}
	catch (const Ort::Exception& e) {
		std::cerr << "error loading the model \n" << e.what() << std::endl;
		return MODEL_LOADING_FAILED;
	}


	return MODEL_LOADED_SUCCESSFULLY;
}

extern "C"
void torch_wrap_predict(float i_buffer[], int i_size, float o_buffer[], int* o_size, int* retval, int loopback) {
	static thread_local bool is_loaded = false;
	bool do_log = true;
	int status = 0;

	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/static/static_model.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/seed_609_grad_reg/export_model_01-04-23.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/seed_345_grad_reg_eq_constraints/export_model_01-20-23.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/seed_345_grad_reg_eq_constraints_v2/export_model_03-21-23.pt";
	const char * script_path = "/g/data/up6/nl4727/torchclim-onnx/torch-wrapper/ml-models/version_13/export_model_11-17-24.onnx";

	//test mode
	if( loopback ) {
		for( int i=0; i < (*o_size); i++) {
			o_buffer[i] = (float)(i+1);
		}

		*retval = 3;
		return;
	}


	if( !is_loaded ) {
		status = torch_wrap_create( script_path );
		if( status != MODEL_LOADED_SUCCESSFULLY ) {
        		(*retval) = -1;
			return;
		}
		is_loaded = true;
	}

	// define shape
	const std::array<int64_t, 2> inputShape = { 1, i_size };
	const std::array<int64_t, 2> outputShape = { 1, *o_size };

	// defining tensors
	auto memory_info = Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU);
	Ort::Value inputTensor = Ort::Value::CreateTensor<float>(memory_info, i_buffer, i_size, inputShape.data(), inputShape.size());
	Ort::Value outputTensor = Ort::Value::CreateTensor<float>(memory_info, o_buffer, *o_size, outputShape.data(), outputShape.size());

	// define names
	Ort::AllocatorWithDefaultOptions ort_alloc;
	Ort::AllocatedStringPtr inputName = session.GetInputNameAllocated(0, ort_alloc);
	Ort::AllocatedStringPtr outputName = session.GetOutputNameAllocated(0, ort_alloc);
	const std::array<const char*, 1> inputNames = { inputName.get() };
	const std::array<const char*, 1> outputNames = { outputName.get() };
	inputName.release();
	outputName.release();

	// Execute the model and turn its output into a tensor.
	try {
		Ort::RunOptions runOptions;
		session.Run(runOptions, inputNames.data(), &inputTensor, 1, outputNames.data(), &outputTensor, 1);
	}
	catch (Ort::Exception& e) {
		std::cerr << e.what() << std::endl;
		return;
	}

	// get the shape info from the output tensor
	Ort::TensorTypeAndShapeInfo tensor_info = outputTensor.GetTensorTypeAndShapeInfo();
	std::vector<int64_t> tensor_shape = tensor_info.GetShape();

        // Calculate the number of elements
        int num_elements = 1;
        for (int dim : tensor_shape) {
            num_elements *= dim;
        }

	if(*o_size != num_elements) {
		std::ostringstream msg;
		msg << "torch-wrap output buffer size (" << *o_size <<
			") does not match the model output size (" << num_elements << ").";
		throw std::invalid_argument(msg.str());
	}

	// not needed, automatically copied over to o_buffer in session.run above
	// std::memcpy(o_buffer, output.data_ptr<float>(), (*o_size)*sizeof(float));


	// function is complete
	(*retval) = 1;

	//look for extreme radiation and precip values
	//for (int i=52; i<68; i++) {
	//	if(output[0][i].item<float>() > 1500.0) {
	//		do_log = true;
	//		break;
	//	}

	//}

	if (do_log) {
		// Prepping values for display
		Ort::TensorTypeAndShapeInfo input_type_info = inputTensor.GetTensorTypeAndShapeInfo();
		int input_total_elements = input_type_info.GetElementCount();
		std::string input_dtype = get_type_str(input_type_info.GetElementType());

		Ort::TensorTypeAndShapeInfo output_type_info = outputTensor.GetTensorTypeAndShapeInfo();
		int output_total_elements = output_type_info.GetElementCount();
		std::string output_dtype = get_type_str(output_type_info.GetElementType());

		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap log start:\n";
		std::cout << "torch-wrap inputs : ";
		print_tensor(inputTensor);
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap outputs: ";
		print_tensor(outputTensor);
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap i_tens num elements: " << input_total_elements << "\n";
		std::cout << "torch-wrap outputs num elements: " << output_total_elements << "\n";
		std::cout << "torch-wrap i_tens elements dtype: " << input_dtype << "\n";
		std::cout << "torch-wrap outputs elements dtype: " << output_dtype << "\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap log end:\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	}
}

// helper functions since there is no neat way to run std::cout << tensor like in pytorch
extern "C"
void print_tensor(Ort::Value& tensor) {
	Ort::TensorTypeAndShapeInfo tensor_info = tensor.GetTensorTypeAndShapeInfo();
	std::vector<int64_t> shape = tensor_info.GetShape();
	std::cout << "Tensor Shape: ";
	int total_size = 1;
	for (const int dim : shape) {
		total_size *= dim;
		std::cout << dim << " ";
	}

	const float* data = tensor.GetTensorData<float>();

	std::cout << "\nTensor Values: ";
	for (int i = 0; i < total_size; ++i) {
		std::cout << data[i] << " ";
	}
	std::cout << std::endl;
}

// helper function as the tensor.GetElementType() does not return a readable string, requires transformation to be readable in the output
extern "C"
std::string get_type_str(ONNXTensorElementDataType dtype) {
	std::string dtype_str;
	switch (dtype) {
	case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT:
		dtype_str = "float";
		break;
	case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT32:
		dtype_str = "int32";
		break;
	case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64:
		dtype_str = "int64";
		break;
	case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT8:
		dtype_str = "uint8";
		break;
	case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT8:
		dtype_str = "int8";
		break;
	case ONNX_TENSOR_ELEMENT_DATA_TYPE_BOOL:
		dtype_str = "bool";
		break;
	case ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE:
		dtype_str = "double";
		break;
	default:
		dtype_str = "unknown";
		break;
	}
	return dtype_str;
}
