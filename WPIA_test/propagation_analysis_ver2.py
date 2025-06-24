import numpy as np
import pytplot
from datetime import datetime
import scipy.signal as signal
import numpy as np
from datetime import datetime, timedelta
import xarray as xr
import numpy as np
import pandas as pd
from pytplot import tplot, data_quants, store_data, options, split_vec, cdf_to_tplot, xlim, get_data
import sys
import os
import numpy as np
import math
import spiceypy as spice
import datetime
from datetime import datetime, timedelta
import pyspedas
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, resample, stft
from scipy.linalg import svd
import pyspedas
from pyspedas import tinterpol

def data_read(sts_file_path, trange=None):
    """ 
    Read the STS file and store the data in pytplot

    Parameters:
    sts_file_path (str): The path to the STS file
    trange (list): The time range to extract data from the STS file
        ex. trange=['2021-06-07 16:30:04', '2021-06-07 17:10:04']
    """
    # Step 1: Read the file and process the lines in chunks for efficiency
    times = []
    magnetic_fields = []

    with open(sts_file_path, 'r') as sts_file:
        lines = sts_file.readlines()[119:]  # Start reading from line 120
        for i, line in enumerate(lines, start=120):
            parts = line.split()
            if i == 150:  # Extract year from line 150
                year = int(parts[0])
            if len(parts) >= 10:  # Only process if the line has enough parts
                times.append(parts[6])
                magnetic_fields.append([float(parts[7]), float(parts[8]), float(parts[9])])

    # Step 2: Convert times to timestamps using vectorized numpy operations
    decimal_days = np.array(times, dtype=float)
    base_date = datetime(int(year), 1, 1)

    # Split the decimal day into whole days and fractional parts
    days = decimal_days.astype(int)
    fractions = decimal_days - days

    # Vectorized timedelta calculation for the fractional day part
    offsets = fractions * 86400  # Convert fractional days into seconds

    # Create timestamps in the desired format including nanoseconds
    timestamps = []
    for day, offset in zip(days, offsets):
        current_time = base_date + timedelta(days=int(day) - 1, seconds=float(offset))  # Convert numpy.int64 to native int
        nanoseconds = int((offset - int(offset)) * 1e9)  # Calculate nanoseconds
        # Format as 'YYYY-MM-DDTHH:MM:SS.ffffff' and append nanoseconds
        timestamp = f'{current_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]}{nanoseconds % 1000:03d}'
        timestamps.append(timestamp)

    # Step 3: Convert timestamps to numpy datetime64 format
    datetime64_timestamps = np.array(timestamps, dtype='datetime64[ns]')

    # Step 4: Store the data using pytplot
    magnetic_fields = np.array(magnetic_fields)  # Convert list to numpy array
    pytplot.store_data('juno_mag_24h', data={'x': datetime64_timestamps, 'y': magnetic_fields})

    juno_mag = data_quants['juno_mag_24h']

    #juno_mag_clean = juno_mag.drop_duplicates(dim='time')

    if trange is not None:
        start_time = np.datetime64(trange[0])
        end_time = np.datetime64(trange[1])

        # 時間の差分を計算し、最も近い値のインデックスを取得
        start_idx = np.abs(juno_mag['time'].values - start_time).argmin()
        end_idx = np.abs(juno_mag['time'].values - end_time).argmin()

        # 最も近いインデックスでデータをスライス
        juno_mag_filtered = juno_mag.isel(time=slice(start_idx, end_idx))

        store_data('juno_mag_pre', data={'x': juno_mag_filtered['time'], 'y': juno_mag_filtered})

    else:
        store_data('juno_mag_pre', data={'x': juno_mag['time'], 'y': juno_mag})
    
    return

def cyclotron_frequency(mass, charge, magnetic_field):
    return charge * magnetic_field / mass / (2 * np.pi)

def calc_cyclotron_freq():
    """ 
    Calculate the cyclotron frequency using the magnetic field data

    Returns:
    float: The cyclotron frequency in Hz
    """

    mag = data_quants['juno_mag_pre'].sortby('time')
    scalars = np.sqrt(mag[:, 0]**2 + mag[:, 1]**2 + mag[:, 2]**2) * 1e-9
    proton_mass = 1.6726219e-27

    oxygen_mass = proton_mass * 16
    h2o_mass = proton_mass * 18
    ho_mass = oxygen_mass + proton_mass
    electron_charge = 1.60217662e-19
    o2_mass = oxygen_mass * 2
    sulfur_mass = proton_mass * 32
    h2 = proton_mass * 2

    proton_cyclotron = cyclotron_frequency(proton_mass, electron_charge, scalars)
    oxygen_cyclotron = cyclotron_frequency(oxygen_mass, electron_charge, scalars)
    ho_cyclotron = cyclotron_frequency(ho_mass, electron_charge, scalars)
    h2o_cyclotron = cyclotron_frequency(h2o_mass, electron_charge, scalars)
    o2_cyclotron = cyclotron_frequency(o2_mass, electron_charge, scalars)
    sulfur_cyclotron = cyclotron_frequency(sulfur_mass, electron_charge, scalars)
    h2_cyclotron = cyclotron_frequency(h2, electron_charge, scalars)

    return proton_cyclotron, oxygen_cyclotron, ho_cyclotron, h2o_cyclotron, o2_cyclotron, sulfur_cyclotron, h2_cyclotron

def calculate_moving_average(mag_origin, new_time, window_length_time):
    half_window = pd.to_timedelta(window_length_time / 2, unit='s')
    averaged_data = []

    for t in new_time.values:
        # 中心を `t` とし、前後 `window_length_time/2` 秒の範囲を設定
        start_time = t - half_window
        end_time = t + half_window

        # 指定範囲のデータを抽出
        window_data = mag_origin.sel(time=slice(start_time, end_time))

        # 範囲内のデータがあれば平均を計算、なければ NaN を設定
        if len(window_data['time']) > 0:
            avg = window_data.mean(dim='time')
            averaged_data.append(avg.values)
        else:
            averaged_data.append(np.full((3,), np.nan))  # v_dim のサイズに合わせて NaN を設定

    # 結果を xarray.DataArray に変換
    averaged_array = xr.DataArray(
        np.array(averaged_data),
        coords={'time': new_time, 'v_dim': mag_origin.coords['v_dim']},
        dims=['time', 'v_dim']
    )

    return averaged_array

def calc_svd(new_sampling_rate, overwrap_rate, window_length):
    """
    Calculate the SVD of the magnetic field data

    Parameters:
    new_sampling_rate (int): The new sampling rate to resample the data
    overwrap_rate (float): The overwrap rate for the STFT
    window_length (int): The window length for the STFT

    Returns:
    please use pytplot.tplot_names() to see the all of tplot variables
    kvec_LASVD: The wave normal angle
    kvecazm_LASVD: The azimuthal wave normal angle
    polarization_LASVD: The polarization
    planarity_LASVD: The planarity
    lambda1_LASVD: The first eigenvalue
    lambda2_LASVD: The second eigenvalue
    lambda3_LASVD: The third eigenvalue
    The analysis results are based on the following paper:
    Santolík, O., Parrot, M., & Lefeuvre, F. (2003). Singular value decomposition methods for wave propagation analysis. Radio Science, 38(1). https://doi.org/10.1029/2000RS002523
    """

    # load the data
    juno_mag = data_quants['juno_mag_pre']

    # resample the data in the new sampling rate
    start_time = np.datetime64(juno_mag.time.values[0])
    end_time = np.datetime64(juno_mag.time.values[-1])

    time_length = end_time - start_time

    time_interval = np.timedelta64(int(1e9 / new_sampling_rate), 'ns')

    num_points = int(time_length / time_interval)

    time_array = np.array([start_time + i * time_interval for i in range(num_points)])

    store_data('juno_mag_time', data={'x': time_array, 'y': time_array})

    tinterpol(names='juno_mag_pre', interp_to='juno_mag_time')

    juno_mag = data_quants['juno_mag_pre-itrp']

    mag_resampled = pytplot.data_quants['juno_mag_pre-itrp']
    time_resampled = mag_resampled['time']
    mag_origin = data_quants['juno_mag_pre']

    moving_avg = calculate_moving_average(mag_origin, time_resampled, window_length*time_interval)


    store_data('mag_interpolated', data={'x': time_resampled, 'y': moving_avg})

    split_vec('mag_interpolated')
    # rotaion matrix
    data_x = pytplot.get_data('mag_interpolated_x')
    data_y = pytplot.get_data('mag_interpolated_y')
    data_z = pytplot.get_data('mag_interpolated_z')

    rotmat=np.zeros((3, 3, (len(data_x[0]))))
    rotmat_t=np.zeros((3, 3, (len(data_x[0]))))

    for i in range(len(data_x[1])-1):
        bvec = [data_x[1][i], data_y[1][i], data_z[1][i]]
        zz = [0, 0, 1]

        yhat = np.cross(zz, bvec)
        xhat = np.cross(yhat, bvec)
        zhat = bvec

        yhat = yhat / np.linalg.norm(yhat)
        xhat = xhat / np.linalg.norm(xhat)
        zhat = zhat / np.linalg.norm(zhat)

        rotmat[:,:,i] = np.array([xhat, yhat, zhat])
        rotmat_t[:,:,i] = np.array([xhat, yhat, zhat]).T

    # rotate moving average data
    for i in range(len(time_resampled)):
        moving_avg[i] = np.dot(rotmat[:,:,i], mag_resampled[i])

    # calc the complex amplitude
    data = moving_avg.values

    data_resampled_x = data[:, 0]
    data_resampled_y = data[:, 1]
    data_resampled_z = data[:, 2]

    frequencies, times, Zxx = stft(data_resampled_x, fs=new_sampling_rate, nperseg=window_length, noverlap=int(window_length*overwrap_rate), window='hann')
    frequencies, times, Zyy = stft(data_resampled_y, fs=new_sampling_rate, nperseg=window_length, noverlap=int(window_length*overwrap_rate), window='hann')
    frequencies, times, Zzz = stft(data_resampled_z, fs=new_sampling_rate, nperseg=window_length, noverlap=int(window_length*overwrap_rate), window='hann')

    time_start = pd.to_datetime(juno_mag.time.values[0])
    ut_times = [time_start + timedelta(seconds=t) for t in times]

    # Calculate conjugates of each component
    Zxx_conj = np.conj(Zxx)
    Zyy_conj = np.conj(Zyy)
    Zzz_conj = np.conj(Zzz)

    # Compute inner products (complex amplitudes and their conjugates)
    # 9 combinations in total
    product_xx_xx = Zxx * Zxx_conj  # Zxx and Zxx conjugate
    product_xx_yy = Zxx * Zyy_conj  # Zxx and Zyy conjugate
    product_xx_zz = Zxx * Zzz_conj  # Zxx and Zzz conjugate

    product_yy_xx = Zyy * Zxx_conj  # Zyy and Zxx conjugate
    product_yy_yy = Zyy * Zyy_conj  # Zyy and Zyy conjugate
    product_yy_zz = Zyy * Zzz_conj  # Zyy and Zzz conjugate

    product_zz_xx = Zzz * Zxx_conj  # Zzz and Zxx conjugate
    product_zz_yy = Zzz * Zyy_conj  # Zzz and Zyy conjugate
    product_zz_zz = Zzz * Zzz_conj  # Zzz and Zzz conjugate

    # Transpose each product so that the first axis becomes time and the second becomes frequency
    product_xx_xx = product_xx_xx.transpose()  # Now time is the first axis, and frequency is the second
    product_xx_yy = product_xx_yy.transpose()
    product_xx_zz = product_xx_zz.transpose()

    product_yy_xx = product_yy_xx.transpose()
    product_yy_yy = product_yy_yy.transpose()
    product_yy_zz = product_yy_zz.transpose()

    product_zz_xx = product_zz_xx.transpose()
    product_zz_yy = product_zz_yy.transpose()
    product_zz_zz = product_zz_zz.transpose()


    real_p = np.real(product_xx_xx)
    imag_p = np.imag(product_xx_xx)
    com_amp_matrix_00 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_xx_yy)
    imag_p = np.imag(product_xx_yy)
    com_amp_matrix_01 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_xx_zz)
    imag_p = np.imag(product_xx_zz)
    com_amp_matrix_02 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_yy_xx)
    imag_p = np.imag(product_yy_xx)
    com_amp_matrix_10 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_yy_yy)
    imag_p = np.imag(product_yy_yy)
    com_amp_matrix_11 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_yy_zz)
    imag_p = np.imag(product_yy_zz)
    com_amp_matrix_12 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_zz_xx)
    imag_p = np.imag(product_zz_xx)
    com_amp_matrix_20 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_zz_yy)
    imag_p = np.imag(product_zz_yy)
    com_amp_matrix_21 = np.stack((real_p, imag_p), axis=-1)

    real_p = np.real(product_zz_zz)
    imag_p = np.imag(product_zz_zz)
    com_amp_matrix_22 = np.stack((real_p, imag_p), axis=-1)


    # dataset for SVD
    rr = np.zeros((3, 3, len(times), len(frequencies), 2))
    rr[0, 0, :, :, :] = com_amp_matrix_00
    rr[0, 1, :, :, :] = com_amp_matrix_01
    rr[0, 2, :, :, :] = com_amp_matrix_02
    rr[1, 0, :, :, :] = com_amp_matrix_10
    rr[1, 1, :, :, :] = com_amp_matrix_11
    rr[1, 2, :, :, :] = com_amp_matrix_12
    rr[2, 0, :, :, :] = com_amp_matrix_20
    rr[2, 1, :, :, :] = com_amp_matrix_21
    rr[2, 2, :, :, :] = com_amp_matrix_22



    store_data('s00', data={'x': ut_times, 'y': ut_times})

    # mag_origin = data_quants['juno_mag_pre']
    # new_time = data_quants['s00']
    # window_length_time = window_length * time_interval
    # new_time = new_time.astype('datetime64[ns]')
    # moving_avg = calculate_moving_average(mag_origin, new_time, window_length_time)
    # store_data('mag_interpolated', data={'x': new_time, 'y': moving_avg})

    # # interpolate 
    # # tinterpol(names='juno_mag_pre', interp_to='s00', newname='mag_interpolated')
    # # split_vec('mag_interpolated')

    # # rotaion matrix
    # data_x = pytplot.get_data('mag_interpolated_x')
    # data_y = pytplot.get_data('mag_interpolated_y')
    # data_z = pytplot.get_data('mag_interpolated_z')

    # rotmat=np.zeros((3, 3, (len(data_x[0]))))
    # rotmat_t=np.zeros((3, 3, (len(data_x[0]))))

    # for i in range(len(data_x[1])-1):
    #     bvec = [data_x[1][i], data_y[1][i], data_z[1][i]]
    #     zz = [0, 0, 1]

    #     yhat = np.cross(zz, bvec)
    #     xhat = np.cross(yhat, bvec)
    #     zhat = bvec

    #     yhat = yhat / np.linalg.norm(yhat)
    #     xhat = xhat / np.linalg.norm(xhat)
    #     zhat = zhat / np.linalg.norm(zhat)

    #     rotmat[:,:,i] = np.array([xhat, yhat, zhat])
    #     rotmat_t[:,:,i] = np.array([xhat, yhat, zhat]).T

    # # rotate the matrix
    # for i in range(len(times)-1):
    #     for j in range(len(frequencies)-1):
    #         for k in range(2):
    #             rr[:,:,i,j,k] = np.dot(rotmat[:,:,i], np.dot(rr[:,:,i,j,k], rotmat_t[:,:,i]))


    # moving average
    if 'moving_average' not in locals():
        moving_average = 3

    if moving_average is None:
        moving_average = 3

    rrr = np.zeros_like(rr)

    for i in range(times.size):
        for j in range(1, frequencies.size-1):
            idx_j = np.arange(j-moving_average//2, j+moving_average//2+1).astype(int) # ある周波数jの周囲の移動平均につかうデータ数分のindexを取得
            idx_j = np.clip(idx_j, 0, frequencies.size-1)  # indexが周波数の取りうる範囲に収まるように

            rrr[0,0,i,j,0] = np.sum(rr[0,0,i,idx_j,0]) / float(moving_average)
            rrr[0,0,i,j,1] = np.sum(rr[0,0,i,idx_j,1]) / float(moving_average)
            rrr[0,1,i,j,0] = np.sum(rr[0,1,i,idx_j,0]) / float(moving_average)
            rrr[0,1,i,j,1] = np.sum(rr[0,1,i,idx_j,1]) / float(moving_average)
            rrr[0,2,i,j,0] = np.sum(rr[0,2,i,idx_j,0]) / float(moving_average)
            rrr[0,2,i,j,1] = np.sum(rr[0,2,i,idx_j,1]) / float(moving_average)
            rrr[1,0,i,j,0] = np.sum(rr[1,0,i,idx_j,0]) / float(moving_average)
            rrr[1,0,i,j,1] = np.sum(rr[1,0,i,idx_j,1]) / float(moving_average)
            rrr[1,1,i,j,0] = np.sum(rr[1,1,i,idx_j,0]) / float(moving_average)
            rrr[1,1,i,j,1] = np.sum(rr[1,1,i,idx_j,1]) / float(moving_average)
            rrr[1,2,i,j,0] = np.sum(rr[1,2,i,idx_j,0]) / float(moving_average)
            rrr[1,2,i,j,1] = np.sum(rr[1,2,i,idx_j,1]) / float(moving_average)
            rrr[2,0,i,j,0] = np.sum(rr[2,0,i,idx_j,0]) / float(moving_average)
            rrr[2,0,i,j,1] = np.sum(rr[2,0,i,idx_j,1]) / float(moving_average)
            rrr[2,1,i,j,0] = np.sum(rr[2,1,i,idx_j,0]) / float(moving_average)
            rrr[2,1,i,j,1] = np.sum(rr[2,1,i,idx_j,1]) / float(moving_average)
            rrr[2,2,i,j,0] = np.sum(rr[2,2,i,idx_j,0]) / float(moving_average)
            rrr[2,2,i,j,1] = np.sum(rr[2,2,i,idx_j,1]) / float(moving_average)

    rr_ = rrr



    n_t = times.size
    n_e = frequencies.size

    A = np.zeros((6, 3, n_t, n_e), dtype=np.float64)
    W2 = np.zeros((3, n_t, n_e), dtype=np.float64)
    V2 = np.zeros((3, 3, n_t, n_e), dtype=np.float64)
    W_SORT = np.zeros((3, n_t, n_e), dtype=np.float64)
    V_SORT = np.zeros((3, 3, n_t, n_e), dtype=np.float64)

    A[0,0,:,:], A[0,1,:,:], A[0,2,:,:] = rr_[0,0,:,:,0], rr_[0,1,:,:,0], rr_[0,2,:,:,0]
    A[1,0,:,:], A[1,1,:,:], A[1,2,:,:] = rr_[1,0,:,:,0], rr_[1,1,:,:,0], rr_[1,2,:,:,0]
    A[2,0,:,:], A[2,1,:,:], A[2,2,:,:] = rr_[2,0,:,:,0], rr_[2,1,:,:,0], rr_[2,2,:,:,0]
    A[3,0,:,:], A[3,1,:,:], A[3,2,:,:] = 0.0, -rr_[0,1,:,:,1], -rr_[0,2,:,:,1]
    A[4,0,:,:], A[4,1,:,:], A[4,2,:,:] = rr_[0,1,:,:,1], 0.0, -rr_[1,2,:,:,1]
    A[5,0,:,:], A[5,1,:,:], A[5,2,:,:] = rr_[0,2,:,:,1], rr_[1,2,:,:,1], 0.0


    A_clean = np.empty_like(A)
    for i in range(n_t):
        for j in range(n_e):
            # Perform SVD on the reshaped slice of A


            """ # データの前処理でNaNやInfを除去
            A_clean[:, :, i, j] = A[:, :, i, j]
            if np.isnan(A_clean[:, :, i, j]).any() or np.isinf(A_clean[:, :, i, j]).any():
                continue  # または適切な処理 """

            A_clean[:, :, i, j] = np.nan_to_num(A[:, :, i, j], nan=0.0, posinf=np.finfo(np.float64).max, neginf=np.finfo(np.float64).min)

            U, W, Vt = svd(A_clean[:, :, i, j], full_matrices=False, lapack_driver='gesvd')
            W2[:, i, j] = W
            V2[:, :, i, j] = Vt#.T[:3, :3]  # Note: Vt.T is the transpose of Vt to match IDL's V

            # Sort singular values and corresponding vectors
            W_order = np.argsort(W2[:, i, j])
            for k in range(3):
                W_SORT[k, i, j] = W2[W_order[k], i, j]
                V_SORT[k, :, i, j] = V2[W_order[k], :, i, j]

    powspec_b = np.zeros((n_t, n_e))
    wna = np.zeros((n_t, n_e))
    wna_azm = np.zeros((n_t, n_e))
    polarization = np.zeros((n_t, n_e))
    planarity = np.zeros((n_t, n_e))
    lambda1 = np.zeros((n_t, n_e))
    lambda2 = np.zeros((n_t, n_e))
    lambda3 = np.zeros((n_t, n_e))

    # Loop over elements
    for i in range(n_t):
        for j in range(n_e):
            # Power spectrum
            powspec_b[i,j] = np.sqrt(A[0,0,i,j]**2 + A[1,1,i,j]**2 + A[2,2,i,j]**2)
            # Wave normal
            ## 0で割る作業に対しては90度を返すように
            if V_SORT[0,2,i,j] == 0.:
                wna[i,j] = 90.
            else:
                wna[i,j] = np.abs(np.arctan(np.sqrt(V_SORT[0,0,i,j]**2+V_SORT[0,1,i,j]**2)/V_SORT[0,2,i,j])/np.pi*180.) #[degree]
            wna_azm[i,j] = np.arctan2(V_SORT[0,1,i,j], V_SORT[0,0,i,j])/np.pi*180.
            # Polarization
            #＃ 0で割る作業に対し値を返さないように設定
            if W_SORT[2,i,j] == 0.:
                polarization[i,j] = np.nan
                planarity[i,j] = np.nan
            else:
                polarization[i,j] = W_SORT[1,i,j]/W_SORT[2,i,j]
                planarity[i,j] = 1. - np.sqrt(W_SORT[0,i,j]/W_SORT[2,i,j])
            if rr_[0, 1, i, j, 1] < 0.:
                polarization[i, j] *= -1.
            # Planarity
            # planarity[i,j] = 1. - np.sqrt(W_SORT[0,i,j]/W_SORT[2,i,j])
            # Lambda
            lambda1[i,j] = W_SORT[0,i,j]
            lambda2[i,j] = W_SORT[1,i,j]
            lambda3[i,j] = W_SORT[2,i,j]

    moving_average = 1
    ma = '{}'.format(new_sampling_rate)
    if moving_average != 1:
        ma = '_ma' + str(moving_average) + '{}'.format(new_sampling_rate)

    pytplot.store_data('powspec_b_LASVD'+ma, data={'x': ut_times, 'y': powspec_b, 'v': frequencies})
    pytplot.options('powspec_b_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'ytitle': 'powspec', 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1, 'ylog':1, 'zlog':1}) 

    pytplot.store_data('kvec_LASVD'+ma, data={'x': ut_times, 'y': wna, 'v': frequencies})
    pytplot.options('kvec_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0.,90.] ,'ytitle': 'WNA', 'ztitle': '[deg]', 'ysubtitle': '[Hz]', 'spec': 1, 'ylog':1})

    pytplot.store_data('kvecazm_LASVD'+ma, data={'x': ut_times, 'y': wna_azm, 'v': frequencies})
    pytplot.options('kvecazm_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0.,360.] ,'ytitle': 'azimuthal WNA', 'ztitle': '[deg]', 'ysubtitle': '[Hz]', 'spec': 1})

    pytplot.store_data('polarization_LASVD'+ma, data={'x': ut_times, 'y': polarization, 'v': frequencies})
    pytplot.options('polarization_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[-1.,1.] ,'ytitle': 'polarization', 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1, 'ylog':1})

    pytplot.store_data('planarity_LASVD'+ma, data={'x': ut_times, 'y': planarity, 'v': frequencies})
    pytplot.options('planarity_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0.,1.] ,'ytitle': 'planarity', 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1, 'ylog':1})

    pytplot.store_data('lambda1_LASVD'+ma, data={'x': ut_times, 'y': lambda1, 'v': frequencies})
    pytplot.options('lambda1_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0,10.] ,'ytitle': 'lambda1!CLA SVD'+ma, 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1})

    pytplot.store_data('lambda2_LASVD'+ma, data={'x': ut_times, 'y': lambda2, 'v': frequencies})
    pytplot.options('lambda2_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0,10.] ,'ytitle': 'lambda2!CLA SVD'+ma, 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1})

    pytplot.store_data('lambda3_LASVD'+ma, data={'x': ut_times, 'y': lambda3, 'v': frequencies})
    pytplot.options('lambda3_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0,10.] ,'ytitle': 'lambda3!CLA SVD'+ma, 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1})

    pytplot.store_data('lambda2-1_LASVD'+ma, data={'x': ut_times, 'y': lambda2-lambda1, 'v': frequencies})
    pytplot.options('lambda2-1_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0,0.1] ,'ytitle': 'lambda2-1!CLA SVD'+ma, 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1})

    pytplot.store_data('lambda3-2_LASVD'+ma, data={'x': ut_times, 'y': lambda3-lambda2, 'v': frequencies})
    pytplot.options('lambda3-2_LASVD'+ma, opt_dict = {'yrange':[0.1, 1], 'zrange':[0,0.1] ,'ytitle': 'lambda3-2!CLA SVD'+ma, 'ztitle': '', 'ysubtitle': '[Hz]', 'spec': 1})

    return

def summary_plot(trange, new_sampling_rate, pwoer_range=None, freq_range=None):
    """
    Plot the summary of the analysis results

    Parameters:
    trange (list): The time range to plot the summary
        ex. trange=['2021-06-07 16:30:04', '2021-06-07 17:10:04']
    new_sampling_rate (int): The new sampling rate to resample the data
    pwoer_range (list): The power range to plot the summary
        ex. pwoer_range=[0.1, 1]
    """
    if pwoer_range is None:
        pwoer_range = [1e-5,1e2]
    if freq_range is None:
        freq_range =[0.03, new_sampling_rate/2]
    xlim(trange[0], trange[1])
    options('powspec_b_LASVD{}'.format(str(new_sampling_rate)), opt_dict = {'yrange':freq_range, 'zrange':pwoer_range})
    options('kvec_LASVD{}'.format(str(new_sampling_rate)), opt_dict = {'yrange':freq_range})
    options('kvecazm_LASVD{}'.format(str(new_sampling_rate)), opt_dict = {'yrange':freq_range})
    options('polarization_LASVD{}'.format(str(new_sampling_rate)), opt_dict = {'yrange':freq_range})
    options('planarity_LASVD{}'.format(str(new_sampling_rate)), opt_dict = {'yrange':freq_range})
    tplot(['powspec_b_LASVD{}'.format(str(new_sampling_rate)), 'kvec_LASVD{}'.format(str(new_sampling_rate)), 'polarization_LASVD{}'.format(str(new_sampling_rate)), 'planarity_LASVD{}'.format(str(new_sampling_rate))], xsize=10, ysize=10)
    return

def main(path):
    # Step 1: Read the STS file and store the data in pytplot
    sts_file_path = path
    trange = ['2021-06-07 16:30:04', '2021-06-07 17:10:04'] # data read time range
    data_read(sts_file_path, trange)

    # Step 2: Calculate the cyclotron frequency using the magnetic field data
    # proton_cyclotron, oxygen_cyclotron, ho_cyclotron, h2o_cyclotron, o2_cyclotron, sulfur_cyclotron, h2_cyclotron = calc_cyclotron_freq()

    # Step 3: Calculate the SVD of the magnetic field data
    new_sampling_rate = 8
    overwrap_rate = 0.93
    window_length = 1024
    calc_svd(new_sampling_rate, overwrap_rate, window_length)

    # Step 4: Plot the summary of the analysis results
    trange = ['2021-06-07 16:30:04', '2021-06-07 17:10:04'] # plot time range
    summary_plot(trange, new_sampling_rate)

    return


sts_file_path = '/mnt/d/JUNO_EMIC/event/fgm_jno_l3_2021158pc_v01.sts'
main(sts_file_path)