import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram


def find_min_distance(matrix):
    min_distance = float('inf')
    min_indices = (0, 0)
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix[i])):
            if matrix[i][j] < min_distance:
                min_distance = matrix[i][j]
                min_indices = (i, j)
    return min_indices, min_distance


def update_distance_matrix(matrix, clusters, pair, method):
    new_matrix = []
    size = len(matrix)

    for i in range(size):
        if i not in pair:
            new_row = []
            for j in range(size):
                if j not in pair:
                    new_row.append(matrix[i][j])
            new_matrix.append(new_row)

    if method == 'Minimum':
        new_row = [min(matrix[pair[0]][k], matrix[pair[1]][k])
                   for k in range(size) if k not in pair]
    elif method == 'Maximum':
        new_row = [max(matrix[pair[0]][k], matrix[pair[1]][k])
                   for k in range(size) if k not in pair]
    elif method == 'Average':
        new_row = [(matrix[pair[0]][k] + matrix[pair[1]][k]) /
                   2 for k in range(size) if k not in pair]

    new_row.append(0.0)
    for row in new_matrix:
        row.append(new_row[new_matrix.index(row)])
    new_matrix.append(new_row)

    return new_matrix


def compute_cophenetic_coefficient(matrix, cophenetic_matrix):
    size = len(matrix)
    sum_matrix, sum_cophenetic = 0.0, 0.0

    for i in range(size):
        for j in range(i + 1, size):
            sum_matrix += matrix[i][j]
            sum_cophenetic += cophenetic_matrix[i][j]

    mean_matrix = sum_matrix / ((size * size - size) / 2)
    mean_cophenetic = sum_cophenetic / ((size * size - size) / 2)

    sum_xy, sum_x_sq, sum_y_sq = 0.0, 0.0, 0.0

    for i in range(size):
        for j in range(i + 1, size):
            x = matrix[i][j] - mean_matrix
            y = cophenetic_matrix[i][j] - mean_cophenetic
            sum_xy += x * y
            sum_x_sq += x ** 2
            sum_y_sq += y ** 2

    if sum_x_sq > 0 and sum_y_sq > 0:
        ccc = sum_xy / (sum_x_sq ** 0.5 * sum_y_sq ** 0.5)
    else:
        ccc = 0.0

    return ccc


def draw_dendrograms(clusters_list, distances_list, methods, save_filename=None):
    plt.figure(figsize=(18, 6))
    for index, (clusters, distances, method) in enumerate(zip(clusters_list, distances_list, methods)):
        Z = []
        cluster_counter = len(clusters) + 1
        for i, (c1, c2) in enumerate(clusters):
            Z.append([c1, c2, distances[i], cluster_counter])
            cluster_counter += 1
        Z = np.array(Z)

        plt.subplot(1, 3, index + 1)
        plt.title(f'Method: {method}')
        plt.xlabel('Cluster')
        dendro = dendrogram(Z, color_threshold=0,
                            link_color_func=lambda k: 'red')

        for i, d, coord in zip(range(len(Z)), Z[:, 2], dendro['icoord']):
            y = d
            plt.text((coord[1] + coord[2]) / 2, y,
                     f'{d:.2f}', ha='center', va='bottom', color='black')

        plt.gca().axes.get_yaxis().set_visible(False)

    if save_filename:
        plt.savefig(save_filename)
    # plt.show()


def hierarchical_clustering(matrix, method):
    clusters = [{i} for i in range(len(matrix))]
    current_matrix = matrix
    result_clusters = []
    result_matrices = [matrix]
    min_distances = []
    cophenetic_matrix = np.zeros((len(matrix), len(matrix)))
    cluster_indices = list(range(len(clusters)))

    while len(clusters) > 1:
        (i, j), min_dist = find_min_distance(current_matrix)
        new_cluster = clusters[i] | clusters[j]

        for m in clusters[i]:
            for n in clusters[j]:
                cophenetic_matrix[m, n] = min_dist
                cophenetic_matrix[n, m] = min_dist

        clusters = [clusters[k]
                    for k in range(len(clusters)) if k != i and k != j]
        clusters.append(new_cluster)
        result_clusters.append((cluster_indices[i], cluster_indices[j]))

        cluster_indices = [cluster_indices[k]
                           for k in range(len(cluster_indices)) if k != i and k != j]
        cluster_indices.append(len(matrix) + len(result_clusters) - 1)

        current_matrix = update_distance_matrix(
            current_matrix, clusters, (i, j), method)
        result_matrices.append(current_matrix.copy())
        min_distances.append(min_dist)

    return result_clusters, result_matrices, min_distances, cophenetic_matrix
