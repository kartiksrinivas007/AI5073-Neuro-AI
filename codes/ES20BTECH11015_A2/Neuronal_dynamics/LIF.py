import numpy as np
import matplotlib.pyplot as plt
import math


class LIF:
    def __init__(self, thresh, tau, E_l, R_l, v_reset, t_ref) -> None:
        """
        Initialize the LIF neuron with the following parameters:
        `tau`: membrane time constant
        `thresh`: threshold
        `E_l`: resting potential
        `R_l`: membrane resistance
        `v_reset`: reset potential
        `t_ref`: refractory period
        """
        self.thresh = thresh
        self.tau = tau
        self.E_l = E_l
        self.R_l = R_l
        self.v_reset = v_reset
        self.t_ref = t_ref
    def simulate(self, T_max, current, dt) -> tuple:
        """
        Run the simulation for a given current and maximum time duration `T_max` with a time step `dt`, returns the voltages, times and spike count
        """
        t = 0
        v_next = 0
        v = 0
        v_0 = 0
        voltages = []
        times = []
        spike_count = 0
        while(t < T_max):
            v_next = v_0 + (self.E_l - v_0 + self.R_l*current)*dt/self.tau
            if(v_next > self.thresh):
                voltages.extend(list(np.zeros(int(self.t_ref/dt))))
                times.extend(list(np.arange(t, t + self.t_ref, dt)))
                v_0 = self.v_reset
                spike_count = spike_count + 1
                t = t + self.t_ref
                continue
            v_0 = v_next
            voltages.append(v_next)
            times.append(t + dt)
            t = t + dt
        return voltages, times, spike_count
    def see_image(self, image, current_response, dt, option="row_wise", time_interval = 1):
        t = 0
        v_next = 0
        v = 0
        v_0 = 0
        voltages = []
        times = []
        spike_count = 0
        T_max = image.reshape(-1).shape[0]* time_interval
        print(T_max)
        if(option == "row_wise"):
            image_reshaped = image.reshape(-1)
        else:
            image_reshaped = image.T.reshape(-1)
        t_last = 0
        count = 0
        while(t < T_max):
            if(t - t_last >= time_interval or t == 0):
                count = count + 1
                t_last = t
            v_next = v_0 + (self.E_l - v_0 + self.R_l*current_response*image_reshaped[count]*dt)/self.tau #(self.E_l - v_0 + self.R_l*current)*dt/self.tau
            if(v_next > self.thresh):
                breakpoint()
                voltages.extend(list(np.zeros(int(self.t_ref/dt))))
                times.extend(list(np.arange(t, t + self.t_ref, dt)))
                v_0 = self.v_reset
                spike_count = spike_count + 1
                t = t + self.t_ref
                continue
            v_0 = v_next
            voltages.append(v_next)
            times.append(t + dt)
            t = t + dt
        return voltages, times, spike_count

class IF:
    def __init__(self, thresh, Cap, E_l, v_reset, t_ref) -> None:
        """
        Initialize the LIF neuron with the following parameters:
        `cap`: membrane capacitance
        `thresh`: threshold
        `E_l`: resting potential
        `v_reset`: reset potential
        `t_ref`: refractory period
        """
        self.thresh = thresh
        self.E_l = E_l
        self.cap = Cap
        self.v_reset = v_reset
        self.t_ref = t_ref
    def simulate(self, T_max, current, dt) -> tuple:
        """
        Run the simulation for a given current and maximum time duration `T_max` with a time step `dt`, returns the voltages, times and spike count
        """
        t = 0
        v_next = 0
        v = 0
        v_0 = 0
        voltages = []
        times = []
        spike_count = 0
        while(t < T_max):
            v_next = v_0 + current*dt/self.cap #(self.E_l - v_0 + self.R_l*current)*dt/self.tau
            if(v_next > self.thresh):
                voltages.extend(list(np.zeros(int(self.t_ref/dt))))
                times.extend(list(np.arange(t, t + self.t_ref, dt)))
                v_0 = self.v_reset
                spike_count = spike_count + 1
                t = t + self.t_ref
                continue
            v_0 = v_next
            voltages.append(v_next)
            times.append(t + dt)
            t = t + dt
        return voltages, times, spike_count
    def see_image(self, image, current_response, dt, option="row_wise", time_interval = 1):
        t = 0
        v_next = 0
        v = 0
        v_0 = 0
        voltages = []
        times = []
        spike_count = 0
        T_max = image.reshape(-1).shape[0]* time_interval
        print(T_max)
        if(option == "row_wise"):
            image_reshaped = image.reshape(-1)
        else:
            image_reshaped = image.T.reshape(-1)
        t_last = 0
        count = 0
        while(t < T_max):
            if(t - t_last >= time_interval or t == 0):
                count = count + 1
                t_last = t
            v_next = v_0 + current_response*image_reshaped[count]*dt/self.cap #(self.E_l - v_0 + self.R_l*current)*dt/self.tau
            if(v_next > self.thresh):
                breakpoint()
                voltages.extend(list(np.zeros(int(self.t_ref/dt))))
                times.extend(list(np.arange(t, t + self.t_ref, dt)))
                v_0 = self.v_reset
                spike_count = spike_count + 1
                t = t + self.t_ref
                continue
            v_0 = v_next
            voltages.append(v_next)
            times.append(t + dt)
            t = t + dt
        return voltages, times, spike_count
    
# def simulate_LIF(T_max:float, current:float, thresh:float, tau:float, E_l:float, R_l:float, v_reset:float, t_ref:float, dt:float):
#     """
#     This function simulates the dynamics of a LIF neuron. using variables like threshold `tau`, `thresh`, `E_l`, `R_l`, `v_reset`, `t_ref` and `dt
#     """
#     t = 0
#     v_next = 0
#     v = 0
#     v_0 = 0
#     voltages = []
#     times = []
#     spike_count = 0
#     while(t < T_max):
#         v_next = v_0 + (E_l - v_0 + R_l*current)*dt/tau
#         if(v_next > thresh):
#             t = t + t_ref
#             v_0 = v_reset
#             spike_count = spike_count + 1
#             continue
#         v_0 = v_next
#         voltages.append(v_next)
#         times.append(t + dt)
#         t = t + dt
#     return voltages, times, spike_count