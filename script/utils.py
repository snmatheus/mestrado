#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pickle
import numpy as np
from constants import path_pos


def load_variable(file_name):
    with open(path_pos + file_name, 'rb') as file:
      return pickle.load(file)


def save_variable(data, file_name):
    with open(path_pos + file_name, 'wb') as file:
      pickle.dump(data, file)


def ln(x):
    from math import log
    return log(x)


def weib(v,k,c):
    return (k/c) * (v/c)**(k-1) * np.exp(-(v/c)**k)


def wpd(ws):
    return 0.5 * 1.225 * ws**3


def extrapolate(ws):
    return ws * ln(100/0.0002)/ln(10/0.0002)


def climatology(datas):
    sum_ = np.zeros(datas['1990'].shape)
    for data in datas.values():
        sum_ += data

    return sum_/len(datas)


def find_nearest(axes, point):
    index = np.abs(axes - point).argmin()
    return index  # If you want show the value -> axes.flat[index]
