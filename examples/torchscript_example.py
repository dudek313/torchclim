import torch

"""The following script demonstrated how to export a PyTorch module into a 
format (torchscript) that can then be loaded by the TorchClim plugin into a 
GCM.

Here we have a "Dummy" torch.nn Module that is converted into torchscript and
saved into a file (dummy.pt). Please refer to the PyTorch / torchscript 
documentation for further details.

"""


class Dummy(torch.nn.Module):
    def __init__(self, n_in,n_out):
        super().__init__()
        self.linear = torch.nn.Linear(n_in,n_out)

    def forward(self,x):
        return self.linear(x)

# create dummy net

n_in = 100
n_out = 50
n_batch = 5

dummy_net = Dummy(100,50)

# create traced script module
dummy_input = torch.randn(n_batch,n_in)
traced_module = torch.jit.trace(dummy_net, dummy_input)
traced_module.save("dummy.pt")

# use dummy module as regular net
device = "cpu" # or cuda:X
dummy_net_torchscript = torch.jit.load("dummy.pt")

assert (dummy_net_torchscript(dummy_input) == dummy_net(dummy_input)).all() #same output
