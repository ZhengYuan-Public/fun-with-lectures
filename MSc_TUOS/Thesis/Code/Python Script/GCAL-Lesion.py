# %%
import numpy as np
import json
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import h5py
from pathlib import Path

base_path = ''


def run_with_parameter(x, y, r):
    def gen_lesion_json(_x, _y, _r, _json_path):
        with open('gcal_base.json', 'r') as f_r:
            data = json.load(f_r)
        data['lesionX'] = _x
        data['lesionY'] = _y
        data['lesionRadius'] = _r

        with open(_json_path, 'w') as f_w:
            json.dump(data, f_w, indent=4)

    folder_path = './data/GCAL-Lesion_' + str(x) + '_' + str(y) + '_' + str(r)
    Path(folder_path).mkdir(parents=True, exist_ok=True)
    file_path = folder_path + '/gcal_lesion_' + str(x) + '_' + str(y) + '_' + str(r) + '.json'
    gen_lesion_json(x, y, r, file_path)

    exec_path = './stevens'
    json_path = file_path
    log_path = folder_path
    para_1 = str(1)
    para_2 = str(1)
    para_3 = str(1)

    sh_cmd = exec_path + ' ' + json_path + ' ' + log_path + ' ' + para_1 + ' ' + para_2 + ' ' + para_3
    os.system(sh_cmd)


def clean_up(_log_path):
    os.chdir(_log_path)
    try:
        Path('images').mkdir(parents=True, exist_ok=True)
        os.system('mv *.png ./images')
    except:
        pass
    try:
        Path('h5_retina').mkdir(parents=True, exist_ok=True)
        os.system('mv retina* ./h5_retina')
    except:
        pass
    try:
        Path('h5_v1').mkdir(parents=True, exist_ok=True)
        os.system('mv v1* ./h5_v1')
    except:
        pass
    try:
        Path('h5_retinotopy').mkdir(parents=True, exist_ok=True)
        os.system('mv retinotopy* ./h5_retinotopy')
    except:
        pass
    try:
        Path('h5_pseudo_stimuli_v1').mkdir(parents=True, exist_ok=True)
        os.system('mv circle_stimulus* ./h5_pseudo_stimuli_v1')
        os.system('mv square_stimulus* ./h5_pseudo_stimuli_v1')
    except:
        pass

    os.chdir(base_path)


def visualize_v1(_h5_path, filename):
    h5_file = h5py.File(_h5_path)
    key_1, key_2, key_3 = h5_file.keys()
    _x = h5_file[f'{key_1}'][:]
    _y = h5_file[f'{key_2}'][:]
    _z = h5_file[f'{key_3}'][:]
    # print("Size of response:", _z.size)
    plt.scatter(_x, _y, c=_z, cmap='jet')
    plt.gca().set_aspect('equal')
    plt.savefig(f'{filename}.eps', format='eps')
    # plt.gca().set_title(f'V1 Response of Pseudo {filename} Stimulus')

    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    plt.show()


def gen_pseudo_retina():
    f_base = h5py.File('pattern_base.h5', 'r')
    _X = f_base['X'][:]
    _Y = f_base['Y'][:]
    f_base.close()

    pseudo_patterns = np.identity(len(_X))

    Path(base_path + '/raw_data/pseudo_retina').mkdir(parents=True, exist_ok=True)
    os.chdir(base_path + '/raw_data/pseudo_retina')

    for i in range(len(_X)):
        f_name = 'pseudo_retina_' + "{:04n}".format(i)
        h5f = h5py.File(f_name, 'w')
        h5f.create_dataset('X', data=_X)
        h5f.create_dataset('Y', data=_Y)
        h5f.create_dataset('v1', data=pseudo_patterns[i])
        h5f.close()

    os.chdir(base_path)


def gen_concentric_circles(circle_samples):
    path_circles = base_path + '/raw_data/pseudo_patterns/circles'
    Path(path_circles).mkdir(parents=True, exist_ok=True)
    os.chdir(path_circles)

    circle_center = (0, 0)

    fig = plt.figure(figsize=(256 / 72, 256 / 72), dpi=72)
    ax = plt.axes()

    for i in range(circle_samples):
        ax.clear()

        radius = i * 0.5 / circle_samples
        circle = plt.Circle(circle_center, radius, color='white', fill=False)
        ax.add_artist(circle)

        low, high = -0.5, 0.5
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.axis("off")
        ax.set_aspect('equal')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        img_name = "{:.3f}".format(radius)
        plt.savefig(f'circle_{img_name}.png', bbox_inches='tight', transparent=True, pad_inches=0)

    os.chdir(base_path)


def gen_concentric_rays(num_samples_per_row):
    path_rays = base_path + '/raw_data/pseudo_patterns/rays'
    Path(path_rays).mkdir(parents=True, exist_ok=True)
    os.chdir(path_rays)

    step_size = 1 / num_samples_per_row
    x_s = [step_size * i for i in range(num_samples_per_row)]
    y_s = [step_size * i for i in range(num_samples_per_row)]

    x_bottom = x_s
    y_bottom = np.zeros(num_samples_per_row)

    x_top = x_s
    y_top = np.zeros(num_samples_per_row) + 1

    x_left = np.zeros(num_samples_per_row)
    y_left = y_s

    x_right = np.zeros(num_samples_per_row) + 1
    y_right = y_s

    X = np.concatenate((x_left, x_right, x_bottom, x_top))
    Y = np.concatenate((y_left, y_right, y_bottom, y_top))

    fig = plt.figure(figsize=(256 / 72, 256 / 72), dpi=72)
    ax = plt.axes()
    ray_center = (0.5, 0.5)

    for index, (_x, _y) in enumerate(zip(X, Y)):
        ax.clear()
        ax.plot((ray_center[0], _x), (ray_center[1], _y), color='white')

        low, high = 0, 1
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.axis("off")
        ax.set_aspect('equal')

        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        img_name = "{:04n}".format(index)
        plt.savefig(f'ray_{img_name}.png', bbox_inches='tight', transparent=True, pad_inches=0)


def gen_shifted_single_pixel_patterns():
    path_rays = base_path + '/raw_data/pseudo_patterns/single_pixel'
    Path(path_rays).mkdir(parents=True, exist_ok=True)
    os.chdir(path_rays)

    img_size = 256
    fig = plt.figure(figsize=(img_size / 72, img_size / 72), dpi=72)
    ax = plt.axes()

    for _x in range(img_size):
        for _y in range(img_size):
            ax.clear()

            low, high = 0, 1
            ax.set_xlim(low, high)
            ax.set_ylim(low, high)
            ax.axis("off")
            ax.set_aspect('equal')

            plt.scatter(_x/img_size, _y/img_size, color='white')

            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            img_name = "{:04n}".format(img_size*_x + _y)
            plt.savefig(f'ray_{img_name}.png', bbox_inches='tight', transparent=True, pad_inches=0)

    os.chdir(base_path)


def visualize_retinotopy(_retinotopy_path):

    retinotopy_file = h5py.File(_retinotopy_path, 'r')
    max_X = retinotopy_file['max_X'][:]
    max_Y = retinotopy_file['max_Y'][:]
    retinotopy_file.close()

    fig = plt.figure()
    ax = plt.axes()
    plt.scatter(max_X, max_Y, s=1)  # s = _markersize_
    low, high = -0.6, 0.6
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_aspect('equal')
    ax.set_title('Retinotopy Map')
    # plt.savefig('retinotopy_map.eps', format='eps')
    plt.show()

    # v1_file = h5py.File(v1_path, 'r')
    # v1_X = v1_file['X'][:]
    # v1_Y = v1_file['Y'][:]
    # v1_response = v1_file['response'][:]
    # print(len(v1_response))
    # v1_file.close()
    # fig2 = plt.figure()
    # ax2 = plt.axes()
    # plt.scatter(v1_X, v1_Y, c=v1_response)
    # low, high = -0.6, 0.6
    # ax2.set_xlim(low, high)
    # ax2.set_ylim(low, high)
    # ax2.set_aspect('equal')
    # ax2.set_title('V1 Response')
    # plt.savefig(f'v1 response@{_time}.eps', format='eps')
    # plt.show()


def visualize_v1_via_retinotopy(_log_path):
    os.chdir(_log_path)

    retinotopy_path = _log_path + '/retinotopy'
    circle_path = _log_path + '/circle'
    square_path = _log_path + '/square'

    retinitopy_h5 = h5py.File(retinotopy_path)
    retinitopy_max_X = retinitopy_h5['max_X'][:]
    retinitopy_max_Y = retinitopy_h5['max_Y'][:]
    retinitopy_h5.close()

    circle_h5 = h5py.File(circle_path)
    circle_X = circle_h5['X'][:]
    circle_Y = circle_h5['Y'][:]
    circle_response = circle_h5['response'][:]
    circle_h5.close()

    square_h5 = h5py.File(square_path)
    square_X = square_h5['X'][:]
    square_Y = square_h5['Y'][:]
    square_response = square_h5['response'][:]
    square_h5.close()

    visualize_v1(circle_path, 'Circle')

    fig = plt.figure()
    ax = plt.axes()
    # low, high = -0.6, 0.6
    # ax.set_xlim(low, high)
    # ax.set_ylim(low, high)
    ax.set_aspect('equal')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.set_title(f'Visualization of Pseudo Circle Stimulus via Retinotopy')

    plt.scatter(retinitopy_max_X, retinitopy_max_Y, s=circle_response * 5)
    plt.savefig('circle_via_retinotopy.eps', format='eps')
    plt.show()

    visualize_v1(square_path, 'Square')

    fig2 = plt.figure()
    ax2 = plt.axes()
    # low, high = -0.6, 0.6
    # ax.set_xlim(low, high)
    # ax.set_ylim(low, high)
    ax2.set_aspect('equal')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.set_title(f'Visualization of Pseudo Square Stimulus via Retinotopy')

    plt.scatter(retinitopy_max_X, retinitopy_max_Y, s=square_response * 2)
    plt.savefig('square_via_retinotopy.eps', format='eps')
    plt.show()

    os.chdir(base_path)



#%%
def random_coordinates(lower_bound, upper_bound):
    r = np.random.uniform(low=lower_bound, high=upper_bound)
    theta = np.random.uniform(low=0, high=360)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    r = np.around(r, 2)
    x = np.around(x, 2)
    y = np.around(y, 2)

    print(r, x, y)

    return r, x, y


r1, x1, y1 = random_coordinates(0, 0.67)
r2, x2, y2 = random_coordinates(0.67, 0.87)
r3, x3, y3 = random_coordinates(0.87, 0.90)
r4, x4, y4 = random_coordinates(0, 0.54)
r5, x5, y5 = random_coordinates(0.54, 0.77)
r6, x6, y6 = random_coordinates(0, 0.42)
r7, x7, y7 = random_coordinates(0.42, 0.65)
