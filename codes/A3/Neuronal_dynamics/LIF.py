import numpy as np
import matplotlib.pyplot as plt


class LIF:
    def __init__(self, thresh, tau, E_l, R_l, v_reset, t_ref ,postsynaptic_neuron) -> None:
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
        self.psn = postsynaptic_neuron
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
                for i in range(int(self.t_ref/dt)):
                    self.psn.update_g_syn("spike", t, v_next)
                # self.psn.times.append(list(np.arange(t, t + self.t_ref, dt)))
                v_0 = self.v_reset
                spike_count = spike_count + 1
                t = t + self.t_ref
                continue
            v_0 = v_next
            voltages.append(v_next)
            times.append(t + dt)
            self.psn.update_g_syn("decay", t, v_next)
            t = t + dt
        return voltages, times, spike_count


class PostSyanpticNeuron:
    def __init__(self, g_synbar,E_syn, tau_syn, dt):
        self.g_synbar = g_synbar
        self.g_syn_t = []
        self.tau_syn = tau_syn
        self.dt = dt
        # self.times = []
        self.synaptic_current = []
        self.E_syn = E_syn
        self.count = 0
    def update_g_syn(self, option, t, voltage=1.0):
        if(t == 0):
            self.g_syn_t.append(0)
            # self.times.append(t)
            self.synaptic_current.append(0)
            self.count = self.count + 1
            # print("hallelujah")
        elif(option == "spike"):
            self.g_syn_t.append(self.g_syn_t[-1] + self.g_synbar*self.dt)
            self.synaptic_current = [(voltage - self.E_syn)*i for i in self.g_syn_t]
            # self.times.append(t)
        elif(option == "decay"):
            self.g_syn_t.append(self.g_syn_t[-1] + -self.g_syn_t[-1]* self.dt/ self.tau_syn)
            if(self.g_syn_t[-1] < 0):
                self.g_syn_t[-1] = 0
                self.count = self.count + 1
            # self.times.append(t)
            self.synaptic_current = [(voltage - self.E_syn)*i for i in self.g_syn_t]
            
class HodgkinHuxley_Neuron():
    """
    Simulate the Hodgkin-Huxley model of a neuron
    """
    def __init__(self, C, g_na, g_k, E_k, E_l,g_l ,E_na):
        self.C = C
        self.g_na = g_na
        self.g_k = g_k
        self.E_k = E_k
        self.E_l = E_l
        self.E_na = E_na
        self.V = -70
        self.g_l  = g_l
        self.alpha_n = lambda v: 0.01*(v + 55)/(1 - np.exp(-(v + 55)/10))
        self.beta_n = lambda v: 0.125*np.exp(-(v + 65)/80)
        self.alpha_m = lambda v: 0.1*(v + 40)/(1 - np.exp(-(v + 40)/10))
        self.beta_m = lambda v: 4*np.exp(-(v + 65)/18)
        self.alpha_h = lambda v: 0.07*np.exp(-(v + 65)/20)
        self.beta_h = lambda v: 1/(1 + np.exp(-(v + 35)/10))
        self.m = 0.05#self.alpha_m(self.V)/(self.alpha_m(self.V) + self.beta_m(self.V))
        self.n = 0.34#self.alpha_n(self.V)/(self.alpha_n(self.V) + self.beta_n(self.V))
        self.h = 0.54#self.alpha_h(self.V)/(self.alpha_h(self.V) + self.beta_h(self.V))

    def step(self, current, dt):
        self.n = self.n + (self.alpha_n(self.V)*(1 - self.n) - self.beta_n(self.V)*self.n)*dt
        self.m = self.m + (self.alpha_m(self.V)*(1 - self.m) - self.beta_m(self.V)*self.m)*dt
        self.h = self.h + (self.alpha_h(self.V)*(1 - self.h) - self.beta_h(self.V)*self.h)*dt
        self.V = self.V + 1/self.C * (current - self.g_na * self.m**3 * self.h * ((self.V + 65) - self.E_na) + self.g_k * self.n**4 * ((self.V + 65) - self.E_k) + self.g_l *((self.V + 65)- self.E_l)) * dt
        return self.V, self.m, self.n, self.h
        pass 
    def simulate(self,T_max, current, dt):
        t = 0
        voltages = []
        ms = []
        ns = []
        hs = []
        times = []
        while (t < T_max):
            v, m , n , h = self.step(current, dt)
            voltages.append(v)
            ms.append(m)
            ns.append(n)
            hs.append(h)
            t = t + dt
            times.append(t)
        pass
        return voltages, ms, ns, hs, times