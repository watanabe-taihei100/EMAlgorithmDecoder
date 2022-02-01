import numpy as np
import math
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas
import subprocess
import random

class DecodeRunner:
    def __init__(self, num_patterns):
        self.num_patterns = num_patterns

        self.Ndx = 896
        self.Ndy = 896
        self.Nmx = 64
        self.Nmy = 64
        self.Nsx = 51
        self.Nsy = 51
        self.mask_pitch = 35
        self.pixel_pitch = 2.5
        self.sky_pitch = 20
        self.detector_to_mask = 25e4

        self.eps = 1e-15
        self.D_minimum = self.eps * np.ones(shape=(self.num_patterns, self.Ndx*self.Ndy), dtype=float)

        self.Mat = None
        self.Pattern = [None]*self.num_patterns
        self.EncodedImage = np.zeros(shape=(self.num_patterns, self.Ndx*self.Ndy))

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

    def construct_matrix(self, attenuation=0.02, narrow=False):
        sky_pos_X = self.detector_to_mask * np.tan(self.arcsec_to_rad(self.sky_X))
        sky_pos_Y = self.detector_to_mask * np.tan(self.arcsec_to_rad(self.sky_Y))

        sky_min_x = np.min(sky_pos_X)
        sky_max_x = np.max(sky_pos_X)
        sky_min_y = np.min(sky_pos_Y)
        sky_max_y = np.max(sky_pos_Y)

        pattern_flatten = [np.append(self.Pattern[pattern_id].flatten().astype(bool), False) for pattern_id in range(self.num_patterns)]

        self.Mat = np.full((self.num_patterns, self.Nsx*self.Nsy, self.Ndx*self.Ndy), attenuation)
        for sky_i, (sky_x, sky_y) in enumerate(zip(sky_pos_X, sky_pos_Y)):
            mask_X = np.round((self.detector_X + sky_x) / self.mask_pitch).astype(int)
            mask_Y = np.round((self.detector_Y + sky_y) / self.mask_pitch).astype(int)

            pattern_index_list = mask_X*self.Nmy+mask_Y
            pattern_index_list[~((mask_X >= 0) & (mask_Y >= 0) & (mask_X < self.Nmx) & (mask_Y < self.Nmy))] = self.Nmx*self.Nmy
            if narrow:
                pattern_index_list[~((self.detector_X >= sky_max_x) & (self.detector_Y >= sky_max_y) & (self.detector_X <= self.pixel_pitch*self.Ndx + sky_min_x) & (self.detector_Y <= self.pixel_pitch*self.Ndy + sky_min_y))] = self.Nmx*self.Nmy
            for pattern_id in range(self.num_patterns):
                self.Mat[pattern_id, sky_i, pattern_flatten[pattern_id][pattern_index_list]] = 1.0
            # normalize
            self.Mat[:, sky_i, :] /= np.sum(self.Mat[:, sky_i, :])
            if sky_i%500 == 0:
                print(sky_i)

    def read_encoded_image_txt(self, filename, pattern_id):
        e_arr = []
        with open(filename, "r") as f:
            for row in f:
                e_arr.append(list(map(int, row.rstrip().split())))
        self.EncodedImage[pattern_id] += np.array(e_arr).flatten()
        print(np.sum(self.EncodedImage))

    def filter_encoded_image(self, number, seed):
        random.seed(seed)
        L = []
        for pattern_id in range(self.num_patterns):
            for ind, num in enumerate(self.EncodedImage[pattern_id]):
                L += [(pattern_id, ind)]*int(num)
        print(len(L), number)
        indexes = random.sample(L, number)
        ret = np.zeros_like(self.EncodedImage)
        for pattern_id, ind in indexes:
            ret[pattern_id, ind] += 1
        self.EncodedImage = ret

    def E_step(self, S_dist):
        return np.dot(S_dist, self.Mat)

    def M_step(self, S_dist, D_dist):
        S_dist_next = np.zeros_like(S_dist, dtype=float)
        for sky_ind in range(self.Nsx*self.Nsy):
            S_dist_next[0, sky_ind] += np.sum(self.EncodedImage * self.Mat[:, sky_ind, :] / np.maximum(D_dist[0], self.D_minimum)) * S_dist[0, sky_ind]

        S_dist_next /= np.sum(S_dist_next)
        return S_dist_next

    def calc_likelihood(self, D_dist):
        L = 0
        for pattern_id in range(self.num_patterns):
            for idx in range(self.Ndx*self.Ndy):
                c = self.EncodedImage[pattern_id, idx]
                if c != 0 and D_dist[0, pattern_id, idx] > 0.0:
                    L += c * D_dist[0, pattern_id, idx]
        return L

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
        return np.loadtxt(filename).T.reshape((1, self.Nsx*self.Nsy))

    def load_divergence_trend(self, filename):
        return list(np.loadtxt(filename))

    def run(self, dirname, max_loop, beta=1, start_index=0):
        # initialize
        if start_index > 0:
            S_dist = self.load_sky_parameters("./parameters/{}/sky_parameters_{}.txt".format(dirname, start_index))
            div_trend = self.load_divergence_trend("./parameters/{}/div_trend.txt".format(dirname))
        else:
            S_dist = np.ones(shape=(1, self.Nsx*self.Nsy), dtype=float) / (self.Nsx*self.Nsy)
            div_trend = []
        D_dist = np.ones(shape=(1, self.num_patterns, self.Ndx*self.Ndy), dtype=float) / (self.num_patterns*self.Ndx*self.Ndy)
        os.makedirs("./pictures/{}".format(dirname), exist_ok=True)
        os.makedirs("./parameters/{}".format(dirname), exist_ok=True)
        print("starting EM algorithm.")
        # loop
        loop_id = start_index
        while loop_id < max_loop:
            print("loop : {} E step".format(loop_id))
            D_dist_t = self.E_step(S_dist)
            print("loop : {} M step".format(loop_id))
            # S_dist = self.M_step(S_dist, D_dist_t)
            S_dist = self.M_step(S_dist, D_dist_t)
            print("finished.")
            likelihood = self.calc_likelihood(D_dist_t)
            div = self.calc_divergence(D_dist, D_dist_t)
            print("likelihood : {}".format(likelihood))
            print("divergence : {}".format(div))
            div_trend.append(div)
            D_dist = D_dist_t

            filename_sky_parameters = "./parameters/{}/sky_parameters_{}.txt".format(dirname, loop_id)
            filename_sky_image = "./pictures/{}/decoded_image_{}.png".format(dirname, loop_id)
            self.save_sky_parameters(filename_sky_parameters, S_dist)
            self.save_divergence_trend("./parameters/{}/div_trend.txt".format(dirname), div_trend)
            # if div < finish_delta:
                # break
            subprocess.run(["../cpp/draw_decoded_image", filename_sky_parameters, filename_sky_image, str(self.Nsx), str(self.Nsy), str(self.sky_pitch)])
            loop_id += 1
        
        print("end.")