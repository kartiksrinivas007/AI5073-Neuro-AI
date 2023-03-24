# %%
import torch
import torchvision
import torchvision.transforms as transforms
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader,Dataset
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
torch.manual_seed(1)
from Neuronal_dynamics import *

# %%
MNIST_data = torchvision.datasets.MNIST(root= './Datasets', download=True, transform = transforms.Compose([transforms.ToTensor(),\
 ]))
data_loader = DataLoader(MNIST_data, batch_size=64, shuffle = True)
test_loader = DataLoader(torchvision.datasets.MNIST('./Datasets', train=False, download=True,
                             transform=torchvision.transforms.Compose([
                               torchvision.transforms.ToTensor(),
                               torchvision.transforms.Normalize(
                                 (0.1307,), (0.3081,))
                             ])), batch_size=10, shuffle = True)
iterator = iter(data_loader)
images, labels = iterator.next()
plt.imshow(images[5].reshape(28, 28), cmap = "gray")
print(images.shape)

# %%
labels

# %%
images, labels = iterator.next()
#convert images to numpy array of size 28 * 28
images  = images.numpy()

# %%
images[11].shape
image = images[11].reshape(28, 28)

# %%
image
image[np.where(image > 0.000)] = 1
image

# %%
plt.imshow(image, cmap = "gray")

# %%
A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
A = A.T.reshape(-1)
print(A)

# %%
E_l = 0
R_l = 1
t_ref = 5
tau = 10 
v_reset = 0
dt = 0.1
T_max = 100
I = 1.2525 # 1.2525 is the best for 0.05Hz freq
v_0 = 0
v_th = 0.96 # 0.96 is the best for 0.05 Hz freq
t = 0
xs = []
ys = []
capacitance = tau/R_l
spike_count = 0
print("capacitance = ", capacitance)

# %%
if_1 = IF(v_th, capacitance, E_l, v_reset, t_ref)
breakpoint()
voltages,times, count = if_1.see_image(image, 100.0, dt, option="row_wise", time_interval = 2)
print("spike_count = ", spike_count)
print(count)

# %%
breakpoint()
len(voltages)
len(times)

# %%
plt.plot(times, voltages)

# %%



