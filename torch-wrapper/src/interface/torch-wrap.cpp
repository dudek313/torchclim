#include <torch/script.h> // One-stop header.

#include <iostream>
#include <memory>
#include <cassert>

static thread_local torch::jit::script::Module the_model;
static thread_local bool is_model_loaded = false;

enum InitState {
	MODEL_LOADING_FAILED = -1,
	MODEL_LOADED_SUCCESSFULLY = 1,
	MODEL_ALREADY_LOADED = 2
};



extern "C"
int torch_wrap_create(const char* script_path) {
	if( is_model_loaded ) {
		return MODEL_ALREADY_LOADED;
	}

	assert( script_path != NULL );

	try {
		// Deserialize the ScriptModule from a file using torch::jit::load().
		the_model = torch::jit::load( script_path );
		the_model.to(torch::kCPU);
		is_model_loaded = true;
	}
	catch (const c10::Error& e) {
		std::cerr << "error loading the model \n" << e.msg();
		return MODEL_LOADING_FAILED;
	}


	return MODEL_LOADED_SUCCESSFULLY;
}

extern "C"
void torch_wrap_predict(float i_buffer[], int i_size, float o_buffer[], int* o_size, int* retval, int loopback) {
	static thread_local torch::jit::script::Module module;
	static thread_local bool is_loaded = false;
	bool do_log = false;

	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch/gaia-deploy-main/traced_model.pt";
	//const char * script_path = "traced_model.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/v2/export_model_cam4_v2.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/v3/export_model_cam4.pt";
	// model v1: const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/v4/export_model_cam4_08-05-22.pt";
	// model v2: 
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/model-2/export_model_08-30-22.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/model-3/export_model_09-01-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/model-4/export_model_09-02-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/v5/export_model_09-15-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/v5/export_model_09-16-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/v5/export_model_09-19-22--no-TS.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_3_100_no_F/export_model_10-12-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_4_100_no_F_no_Z3/export_model_10-12-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_5_100_no_F_no_Z3_positive_PREC/export_model_10-17-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_6_100_no_F_positive_PREC_rect/export_model_10-18-22.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_7_100_no_F_no_Z3_positive_PREC_rect/export_model_10-18-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_8_100_no_F_positive_PREC_softplus/export_model_10-18-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_9_100_no_F_all_positive_softplus/export_model_10-18-22.pt";
	//--const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_9b_100_no_F_all_positive_rect/export_model_10-19-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_9b_100_no_F_all_positive_rect_seed123/export_model_11-21-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_9b_100_no_F_all_positive_rect_seed123/export_seed_345.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_9b_100_no_F_all_positive_rect_seed123/export_seed_609.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_9b_100_no_F_all_positive_rect_seed123/export_seed_112.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_9b_100_no_F_all_positive_rect_1000eps/export_model_11-21-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_A_RELHUM_Q//export_model_11-03-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_B/export_model_11-11-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/version_B_relhum_reg/export_model_11-16-22.pt";
	//const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/static/static_model.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/seed_609_grad_reg/export_model_01-04-23.pt";
	//-const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/seed_345_grad_reg_eq_constraints/export_model_01-20-23.pt";
	const char * script_path = "/g/data/w42/daf561/repo/spcam-ml/libtorch-plugin/torch-wrapper/data/seed_345_grad_reg_eq_constraints_v2/export_model_03-21-23.pt";

	//

	//test mode
	if( loopback ) {
		for( int i=0; i < (*o_size); i++) {
			o_buffer[i] = (float)(i+1);
		}

		*retval = 3;
		return;
	}


	if( !is_loaded ) {
		try {
			// Deserialize the ScriptModule from a file using torch::jit::load().
			module = torch::jit::load( script_path );
			module.to(torch::kCPU);
			is_loaded = true;
		}
		catch (const c10::Error& e) {
		std::cerr << "error loading the model 111\n" << e.msg();
        		(*retval) = -1;
		return;
  		}
	}


	auto options = torch::TensorOptions().dtype(at::kFloat).device(torch::kCPU);
	torch::Tensor i_tens = torch::from_blob(i_buffer, {1, i_size}, options).to(torch::kFloat32);

	std::vector<torch::jit::IValue> inputs;
	/*inputs.push_back(torch::ones({20, 164}));*/
	inputs.push_back( i_tens );

	// Execute the model and turn its output into a tensor.
	at::Tensor output = module.forward(inputs).toTensor();

	if(!output.is_contiguous()) {
		//moving to contigeous memory before copy
		output = output.contiguous();
	}

	//check if Tensor::numel() is a better option - especially if moving to vectors
	if(*o_size != output.sizes()[1]) {
		std::ostringstream msg;
		msg << "torch-wrap output buffer size (" << *o_size << 
			") does not match the model output size (" << output.sizes()[1] << ").";

		output.dtype();
		throw std::invalid_argument(msg.str());
	}

	std::memcpy(o_buffer, output.data_ptr<float>(), (*o_size)*sizeof(float));

	(*retval) = 1;

	//look for extreme radiation and precip values
	//for (int i=52; i<68; i++) { 
	//	if(output[0][i].item<float>() > 1500.0) {
	//		do_log = true;
	//		break;
	//	}
	//}

	if(do_log) {
		//*
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap log start:\n";
		//printf ("psl, srfrad: %6.10f, %6.10f \n", o_buffer[219], o_buffer[214]);
		//std::cout << "torch-wrap output.is_contiguous : " << output.is_contiguous() << "\n";
		//std::cout << "torch-wrap output.is_contiguous : " << output.is_contiguous() << "\n";
		//std::cout << "torch-wrap output.strides : " << output.strides() << "\n";
		//std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap inputs : " << i_tens << "\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap outputs: " << output << "\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap i_tens num elements: " << i_tens.numel() << "\n";
		std::cout << "torch-wrap outputs num elements: " << output.numel() << "\n";
		std::cout << "torch-wrap i_tens elements dtype: " << torch::typeMetaToScalarType(i_tens.dtype()) << "\n";
		std::cout << "torch-wrap outputs elements dtype: " << torch::typeMetaToScalarType(output.dtype()) << "\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "torch-wrap log end:\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		//*/
	}
}


