import numpy as np
import math
import random
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas
import subprocess

class DecodeRunnerWithPolarization:
    def __init__(self, num_patterns):
        self.num_patterns = num_patterns

        self.Ndx = 896
        self.Ndy = 896
        self.Ndp = 2
        self.Nmx = 64
        self.Nmy = 64
        self.Nsx = 51
        self.Nsy = 51
        self.mask_pitch = 35
        self.pixel_pitch = 2.5
        self.sky_pitch = 20
        self.detector_to_mask = 25e4

        self.modulation_factor = 0.10329176923032725
        # self.modulation_factor = 1.0

        self.angle_pitch = 180 / self.Ndp
        self.eps = 1e-15
        self.D_minimum = np.full((self.num_patterns, self.Ndx*self.Ndy), self.eps, dtype=float)

        self.Mat = None
        self.Pattern = [None]*self.num_patterns
        self.EncodedImage = np.zeros(shape=(self.num_patterns, self.Ndp, self.Ndx*self.Ndy))

        self.detector_X = np.array([i//self.Ndy * self.pixel_pitch for i in range(self.Ndy*self.Ndx)])
        self.detector_Y = np.array([i%self.Ndy * self.pixel_pitch for i in range(self.Ndy*self.Ndx)])
        self.sky_X = np.array([(i//self.Nsy - self.Nsx//2) * self.sky_pitch for i in range(self.Nsy*self.Nsx)])
        self.sky_Y = np.array([(i%self.Nsy - self.Nsy//2) * self.sky_pitch for i in range(self.Nsy*self.Nsx)])

    def read_pattern_file(self, filename, pattern_id):
        if not 0 <= pattern_id < self.num_patterns:
            raise Exception("invalid pattern id")
        mask = []
        with open(filename, "r") as f:
            for row in f:
                L = list(map(int, row.rstrip().split()))
                mask.append(L)
        rev_mask = mask[::-1]
        self.Pattern[pattern_id] = np.array(rev_mask)

    def arcsec_to_rad(self, arcsec):
        return arcsec / 3600 / 180 * math.pi

    def construct_matrix(self, attenuation=0.02):
        print("starting construct matrix..")
        sky_pos_X = self.detector_to_mask * np.tan(self.arcsec_to_rad(self.sky_X))
        sky_pos_Y = self.detector_to_mask * np.tan(self.arcsec_to_rad(self.sky_Y))

        pattern_flatten = [np.append(self.Pattern[pattern_id].flatten().astype(bool), False) for pattern_id in range(self.num_patterns)]

        self.Mat = np.full((self.num_patterns, self.Nsx*self.Nsy, self.Ndx*self.Ndy), attenuation)
        for sky_i, (sky_x, sky_y) in enumerate(zip(sky_pos_X, sky_pos_Y)):
            mask_X = np.round((self.detector_X + sky_x) / self.mask_pitch).astype(int)
            mask_Y = np.round((self.detector_Y + sky_y) / self.mask_pitch).astype(int)

            pattern_index_list = mask_X*self.Nmy+mask_Y
            pattern_index_list[~((mask_X >= 0) & (mask_Y >= 0) & (mask_X < self.Nmx) & (mask_Y < self.Nmy))] = self.Nmx*self.Nmy
            for pattern_id in range(self.num_patterns):
                self.Mat[pattern_id, sky_i, pattern_flatten[pattern_id][pattern_index_list]] = 1.0
            # normalize
            self.Mat[:, sky_i, :] /= np.sum(self.Mat[:, sky_i, :])
            if sky_i%500 == 0:
                print(sky_i)
    
    def sampling_data(self, arr, number, seed):
        random.seed(seed)
        L = []
        for ind, num in enumerate(arr):
            L += [ind]*num
        print(len(L), number)
        indexes = random.sample(L, number)
        ret = np.zeros_like(arr)
        for ind in indexes:
            ret[ind] += 1
        return ret

    def read_encoded_image_txt(self, filename, pattern_id, angle, weight=1.0, sample=0, seed=0):
        angle_ind = int(angle / self.angle_pitch + 0.5) % self.Ndp
        e_arr = []
        with open(filename, "r") as f:
            for row in f:
                e_arr.append(list(map(int, row.rstrip().split())))
        arr = np.array(e_arr).flatten()
        if sample > 0:
            arr = self.sampling_data(arr, sample, seed)
        self.EncodedImage[pattern_id, angle_ind] += arr / weight

    def read_encoded_image_txt_double(self, filenamebase, pattern_id, weight_Htype=1.0, weight_Vtype=1.0, sample=0, seed=0):
        e_arr = []
        with open(filenamebase + "_Htype/pattern_{}.txt".format(pattern_id), "r") as f:
            for row in f:
                e_arr.append(list(map(int, row.rstrip().split())))
        arr_Htype = np.array(e_arr).flatten()
        e_arr = []
        with open(filenamebase + "_Vtype/pattern_{}.txt".format(pattern_id), "r") as f:
            for row in f:
                e_arr.append(list(map(int, row.rstrip().split())))
        arr_Vtype = np.array(e_arr).flatten()

        self.EncodedImage[pattern_id, 0] += arr_Htype / weight_Htype
        self.EncodedImage[pattern_id, 1] += arr_Vtype / weight_Vtype
        print(filenamebase)
        print(np.sum(arr_Htype), np.sum(arr_Vtype))


    def E_step(self, S_dist_0, S_dist_90):
        D_dist = np.zeros(shape=(self.num_patterns, self.Ndp, self.Ndx*self.Ndy), dtype=float)
        for sky_i, (s_0, s_90) in enumerate(zip(S_dist_0, S_dist_90)):
            d_Htype = ((1 + self.modulation_factor) * s_0 + (1 - self.modulation_factor) * s_90) / self.Ndp
            d_Vtype = ((1 - self.modulation_factor) * s_0 + (1 + self.modulation_factor) * s_90) / self.Ndp
            for pattern_id in range(self.num_patterns):
                D_dist[pattern_id, 0, :] += self.Mat[pattern_id, sky_i, :] * d_Htype
                D_dist[pattern_id, 1, :] += self.Mat[pattern_id, sky_i, :] * d_Vtype
        norm = np.sum(D_dist)
        print(norm)
        D_dist /= norm
        return D_dist

    def M_step(self, S_dist_0, S_dist_90, D_dist):
        S_dist_0_next = np.zeros_like(S_dist_0, dtype=float)
        S_dist_90_next = np.zeros_like(S_dist_90, dtype=float)

        for sky_i, (s_0, s_90) in enumerate(zip(S_dist_0, S_dist_90)):
            S_dist_0_next[sky_i] += np.sum(self.EncodedImage[:, 0, :] * self.Mat[:, sky_i, :] / np.maximum(D_dist[:, 0, :], self.D_minimum)) * s_0
            S_dist_90_next[sky_i] += np.sum(self.EncodedImage[:, 1, :] * self.Mat[:, sky_i, :] / np.maximum(D_dist[:, 1, :], self.D_minimum)) * s_90
        norm = np.sum(S_dist_90_next) + np.sum(S_dist_0_next)
        S_dist_0_next /= norm
        S_dist_90_next /= norm
        return S_dist_0_next, S_dist_90_next

    def calc_likelihood(self, D_dist):
        return np.sum(self.EncodedImage * D_dist)

    def calc_divergence(self, D_dist_0, D_dist_1):
        D = D_dist_1 - D_dist_0
        return np.sqrt(np.sum(D*D))


    def save_sky_img(self, filename, S_dist):
        X = sorted(list(set(self.sky_X)))
        Y = sorted(list(set(self.sky_Y)))
        df = pandas.DataFrame(np.flip(S_dist.reshape((self.Nsx, self.Nsy)).T), columns=[str(int(x)) for x in X], index=[str(int(y)) for y in reversed(Y)])

        fig = plt.figure(filename)
        sns.heatmap(df, cmap="YlOrRd")
        plt.savefig(filename)
        plt.close(fig)
    
    def save_detector_image(self, filename, D_dist):
        X = sorted(list(set(self.detector_X)))
        Y = sorted(list(set(self.detector_X)))
        df = pandas.DataFrame(np.flip(D_dist.reshape((self.Ndx, self.Ndy)).T), columns=[str(int(x)) for x in X], index=[str(int(y)) for y in reversed(Y)])

        fig = plt.figure(filename)
        sns.heatmap(df, cmap="YlOrRd")
        plt.savefig(filename)
        plt.close(fig)

    def save_sky_parameters(self, filename, S_dist):
        np.savetxt(filename, S_dist.reshape((self.Nsx, self.Nsy)).T)

    def save_detector_parameters(self, filename, D_dist):
        np.savetxt(filename, D_dist.reshape((self.Ndx, self.Ndy)).T)
    
    def save_divergence_trend(self, filename, div_trend):
        np.savetxt(filename, np.array(div_trend))

    def load_sky_parameters(self, filename):
        return np.loadtxt(filename).T.reshape((self.Nsx*self.Nsy, ))

    def load_divergence_trend(self, filename):
        return list(np.loadtxt(filename))

    def run(self, dirname, max_loop, start_index=0):
        print("initialize")
        # initialize
        if start_index == 0:
            S_dist_0 = np.full((self.Nsx*self.Nsy,), 0.5 / (self.Nsx*self.Nsy), dtype=float) 
            S_dist_90 = np.full((self.Nsx*self.Nsy,), 0.5 / (self.Nsx*self.Nsy), dtype=float)
            div_trend = []
        else:
            S_dist_0 = self.load_sky_parameters("./parameters/{}/sky_parameters_{}_H.txt".format(dirname, start_index-1))
            S_dist_90 = self.load_sky_parameters("./parameters/{}/sky_parameters_{}_V.txt".format(dirname, start_index-1))
            div_trend = self.load_divergence_trend("./parameters/{}/div_trend.txt".format(dirname))

        D_dist = np.full((self.num_patterns, self.Ndp, self.Ndx*self.Ndy), 1 / (self.num_patterns*self.Ndx*self.Ndy*self.Ndp), dtype=float)

        os.makedirs("./pictures/{}".format(dirname), exist_ok=True)
        os.makedirs("./parameters/{}".format(dirname), exist_ok=True)
        print("starting EM algorithm.")
        # loop
        loop_id = start_index
        while loop_id < max_loop:
            print("loop : {} E step".format(loop_id))
            D_dist_t = self.E_step(S_dist_0, S_dist_90)
            print("loop : {} M step".format(loop_id))
            S_dist_0, S_dist_90 = self.M_step(S_dist_0, S_dist_90, D_dist_t)
            print("finished.")
            likelihood = self.calc_likelihood(D_dist_t)
            div = self.calc_divergence(D_dist, D_dist_t)
            print("likelihood : {}".format(likelihood))
            print("divergence : {}".format(div))
            div_trend.append(div)
            D_dist = D_dist_t

            filename_Htype = "./parameters/{}/sky_parameters_{}_H.txt".format(dirname, loop_id)
            filename_Vtype = "./parameters/{}/sky_parameters_{}_V.txt".format(dirname, loop_id)
            outfilename_I = "./pictures/{}/decoded_image_{}.png".format(dirname, loop_id)
            outfilename_P = "./pictures/{}/decoded_image_pol_{}.png".format(dirname, loop_id)
            self.save_sky_parameters(filename_Htype, S_dist_0)
            self.save_sky_parameters(filename_Vtype, S_dist_90)
            self.save_divergence_trend("./parameters/{}/div_trend.txt".format(dirname), div_trend)

            subprocess.run(["/Users/watanabe/work/Coded_aperture/cpp/draw_decoded_image_polarization", filename_Htype, filename_Vtype, outfilename_I, outfilename_P, str(self.Nsx), str(self.Nsy), str(self.sky_pitch), str(1.0)])
            loop_id += 1
        
        print("end.")