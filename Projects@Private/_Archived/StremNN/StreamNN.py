import random
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats
import scipy.signal


# -------------Extract Data-------------


def extract_mnist_data(num_1, num_2):
    from mlxtend.data import loadlocal_mnist
    images, labels = loadlocal_mnist(
        images_path='./data/train-images-idx3-ubyte',
        labels_path='./data/train-labels-idx1-ubyte')

    for num in range(10):
        var_name = 'images_' + str(num) + '_index'
        locals()[var_name] = []

    for label_index in range(len(labels)):
        locals()['images_' + str(labels[label_index]) + '_index'].append(label_index)

    index_num_1 = locals()['images_' + str(num_1) + '_index']
    index_num_2 = locals()['images_' + str(num_2) + '_index']

    return images[index_num_1, :], images[index_num_2, :]


# -------------Utilises-------------


def plot_list_list_tuple(_llt):
    def reconstruct_llt(llt_type_data):
        reconstructed_llt = []
        for i in range(len(llt_type_data)):
            for _tuple in llt_type_data[i]:
                reconstructed_llt.append(_tuple)
        return reconstructed_llt
    plot_data = reconstruct_llt(_llt)
    _x, _y = zip(*plot_data)
    plt.plot(_x, _y)
    plt.show()


def plot_list_tuple(_lt):
    plot_data = _lt
    x, y = zip(*plot_data)
    plt.plot(x, y)
    plt.show()


# -------------Neuron Functions-------------


def sequences_in_train(input_intensity, time_in_secs):
    """
    Sequences generated in given training time t(s).
    Output:
    [
        [()()()...(time, 0 or 1) in 1st sec],
        [()()()...(time, 0 or 1) in 2nd sec],
        ...
        [()()()...(time, 0 or 1) in <x> sec],
    ]
    """
    def seq_in_1s(input_intensity_i, time_offset):
        # Sequence generated in 1 sec.
        threshold = 0.8  # f(input_intensity)

        firing_rate_max = 200
        firing_rate_min = 1
        range_fr = firing_rate_max - firing_rate_min

        intensity_max = 255
        intensity_min = 0
        range_intensity = intensity_max - intensity_min

        power = 2

        firing_rate = range_fr * pow((input_intensity_i / range_intensity), power)  # fr = a*(bI)^n

        samples_num = math.ceil(firing_rate) * 2

        seq_generated = []
        for sample in range(samples_num):
            temp = random.random()
            if temp >= threshold:
                seq_generated.append(1)
            else:
                seq_generated.append(0)
        time_axis = np.linspace(0, 1, samples_num) + time_offset
        # plt.scatter(time_axis, seq_generated)
        # plt.show()
        return list(zip(time_axis, seq_generated))
    seqs_list = []
    for i in range(time_in_secs):
        seqs_list.append(seq_in_1s(input_intensity, i))
    return seqs_list


def hot_times_intervals_batches(all_sequence_in_train):
    # Batch = (hot_time, interval between 2 hot_times)
    def calculate_hot_times(sequence_in_train):
        timing = []
        for seq_i in range(len(sequence_in_train)):
            for seq_j in range(len(sequence_in_train[seq_i])):
                if sequence_in_train[seq_i][seq_j][1] == 1:
                    timing.append(sequence_in_train[seq_i][seq_j])
        return timing

    def calculate_intervals(hot_times_sequence):
        intervals_list = []
        for interval_i in range(len(hot_times_sequence)):
            if interval_i == len(hot_times_sequence) - 1:
                t_last = hot_times_sequence[len(hot_times_sequence) - 1][0]
                t_end = math.ceil(t_last)
                intervals_list.append(t_end - t_last)
            else:
                time_1 = hot_times_sequence[interval_i][0]
                time_2 = hot_times_sequence[interval_i + 1][0]
                intervals_list.append(time_2 - time_1)
        return intervals_list

    hot_times_list = calculate_hot_times(all_sequence_in_train)
    interval_list = calculate_intervals(hot_times_list)

    hot_timing = []
    for i in range(len(hot_times_list)):
        hot_timing.append(hot_times_list[i][0])
    return list(zip(hot_timing, interval_list))


def calculate_seqs_effect_receptor(hot_time_and_interval_batches):
    def cal_weight_increment(increment_duration):
        def inc_func(inc_var):
            return pow(inc_var, 2)

        increment = inc_func(increment_duration)
        return increment

    def cal_weight_decrement(decrement_duration):
        def dec_func(dec_var):
            return math.log(dec_var + 1)

        decrement = dec_func(decrement_duration)
        return decrement

    weight_history = []
    initial_weight = random.random()
    current_weight = initial_weight

    t_inc_max = 1 / 100
    max_increment = cal_weight_increment(t_inc_max)

    for i in range(len(hot_time_and_interval_batches)):
        # Weight is updated after each interval.
        batch_i = hot_time_and_interval_batches[i]
        interval_i = batch_i[1]

        if interval_i <= t_inc_max:
            inc_interval = batch_i[1]
            w_inc = cal_weight_increment(inc_interval)
            current_weight += w_inc
            if current_weight >= 1:
                current_weight = 1
        else:
            dec_interval = interval_i - t_inc_max
            w_dec = cal_weight_decrement(dec_interval)
            current_weight += max_increment - w_dec
            if current_weight <= 0:
                current_weight = 0

        weight_history.append(current_weight)
    out_seq_times, _ = zip(*hot_time_and_interval_batches)
    output_sequence = list(zip(out_seq_times, weight_history))
    # weight_history.pop(0)
    return current_weight, weight_history, output_sequence


def weighted_seqs_output():
    pass


def receptor_processor(single_numbers_in_list):
    # Input: [ 0, 244, 224, 212, ...]
    sequences_output_llt = []
    # [ [seq1] [seq2] [seq3] ... [] [] [] [] ]
    for i in range(len(single_numbers_in_list)):
        pixel_i = single_numbers_in_list[i]
        sequence_pixel_i = sequences_in_train(pixel_i, total_train_time)
        batches_pixel_i = hot_times_intervals_batches(sequence_pixel_i)
        _cw, _wh, _out_seq = calculate_seqs_effect_receptor(batches_pixel_i)
        # Current weight, Weight history, Output sequence
        sequences_output_llt.append(_out_seq)
    return sequences_output_llt


def map_layer_i_to_next(_layer_size, _kernel_size, _input_llt):
    """
    Input data:
    [
    [(x1, y1), (x2,y2), (x3,y3),...],
    ... ,
    [(x1, y1), (x2,y2), (x3,y3),...]
    ]
    Output data:
    [
    [[(x1, y1), (x2,y2), (x3,y3),...], [(x1, y1), (x2,y2), (x3,y3),...], ...],
    [[(x1, y1), (x2,y2), (x3,y3),...], [(x1, y1), (x2,y2), (x3,y3),...], ...],
    ...
    ]
    """
    def generate_block(block_start_index):
        indices = []
        for _i in range(_kernel_size):
            for _j in range(_kernel_size):
                index_ij = _layer_size * _i + block_start_index + _j
                indices.append(index_ij)
        block = []
        for block_index in indices:
            block.append(_input_llt[block_index])
        return block

    # Calculate block start index(Upper left).
    rows = int(_layer_size / _kernel_size)
    columns = int(_layer_size / _kernel_size)
    index_0s = []
    for row in range(rows):
        for column in range(columns):
            index_0s.append(row * _layer_size * _kernel_size + column * _kernel_size)

    # Generate all blocks
    blocks = []
    for index_0 in index_0s:
        block_i = generate_block(index_0)
        print(block_i)
        blocks.append(block_i)
    return blocks


def block_summation(out_seqs_from_previous_layer):
    # Input is sequences for post-synaptic neuron.
    # [
    #   [ seq_1],
    #   [...],
    #   [seq_n],
    # ]
    # Output is a single combined sequence. [ (), (), (), ()...]
    temp_sum = []
    for i in range(len(out_seqs_from_previous_layer)):
        for j in range(len(out_seqs_from_previous_layer[i])):
            temp_sum.append(out_seqs_from_previous_layer[i][j])
    temp_sum.sort(key=lambda tup: tup[0])

    def combine_same_t(tp_1, tp_2):
        new_value = tp_1[1] + tp_2[1]
        tp = (tp_1[0], new_value)
        return tp

    for k in range(len(temp_sum)-1, -1, -1):
        if temp_sum[k][0] == temp_sum[k-1][0]:
            _add = combine_same_t(temp_sum[k], temp_sum[k-1])
            temp_sum[k] = _add
            temp_sum.pop(k-1)
    return temp_sum


def activation_function(input_sequences):
    return [(0, 1), (1, 3)]


def neuron_processor(map_ij):
    """
    Input:
    [ [block1], [block2], [block3], ...]
    """
    layer_j_input_seq = []
    # [ [block_sum1], [block_sum2], [block_sum3], ... ]
    for _map in map_ij:
        layer_j_input_seq.append(block_summation(_map))

    layer_j_output_seq = []
    # [ [seq1*weight], [seq2], [seq3], ... ]
    for _seq in layer_j_input_seq:
        layer_j_output_seq.append(activation_function(_seq))

    return layer_j_output_seq


# ------------- Neuron Main -------------


# Extract Data
images_of_2, images_of_4 = extract_mnist_data(2, 4)

# Visualize Data
raw_input = images_of_2[0, :]
raw = raw_input.reshape(28, 28)
plt.imshow(raw, cmap='gray')
plt.show()
plt.savefig('image.png')

# Super Parameters
total_train_time = 50
neurons_in_layer = [28, 14, 7, 1]
kernels_size = [2, 2, 7]

# Start
layer_0_out_seqs = receptor_processor(raw_input)  # [ [seq1] [seq2] [seq3] ... [] [] [] [] ]
map_01 = map_layer_i_to_next(neurons_in_layer[0], kernels_size[0], layer_0_out_seqs)  # [ [Block1] [Block2] ... ]


# ------------- DEBUG -------------

# dummy_data = [
#     20, 30, 234, 21,
#     20, 30, 234, 21,
#     20, 30, 234, 21,
#     20, 30, 234, 21,
# ]
#
# # [
# #     20, 30, 234, 21, 42, 12,
# #     20, 30, 234, 21, 42, 12,
# #     20, 30, 234, 21, 42, 12,
# #     20, 30, 234, 21, 42, 12,
# #     20, 30, 234, 21, 42, 12,
# #     20, 30, 234, 21, 42, 12,
# # ]
#
# dummy_seq_layer0 = neuron_processor(dummy_data)
# layer_1_input = map_layer_i_to_next(4, 4, dummy_seq_layer0)
# temp_in = []
# for i in range(len(layer_1_input)):
#     _sum = synapses_summation(layer_1_input[i])
#     temp_in.append(_sum)