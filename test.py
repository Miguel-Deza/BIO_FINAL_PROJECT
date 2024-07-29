import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

dist_matrix = np.array([
    [0.0, 10.387561, 17.558267, 19.330544, 8.96109, 28.789788, 26.874964, 50.900427, 11.312785, 5.89625],
    [10.387561, 0.0, 20.305233, 24.120915, 13.176571, 35.350279, 27.862574, 56.444247, 12.203479, 11.005102],
    [17.558267, 20.305233, 0.0, 17.94159, 16.932385, 40.157842, 39.459878, 62.46387, 14.361582, 17.345279],
    [19.330544, 24.120915, 17.94159, 0.0, 23.839085, 42.807043, 35.569798, 64.637755, 17.071826, 22.282908],
    [8.96109, 13.176571, 16.932385, 23.839085, 0.0, 26.501404, 29.805982, 49.036852, 12.94894, 6.36347],
    [28.789788, 35.350279, 40.157842, 42.807043, 26.501404, 0.0, 32.249749, 23.603647, 35.288888, 28.741118],
    [26.874964, 27.862574, 39.459878, 35.569798, 29.805982, 32.249749, 0.0, 46.123421, 26.24005, 29.904652],
    [50.900427, 56.444247, 62.46387, 64.637755, 49.036852, 23.603647, 46.123421, 0.0, 56.892871, 51.25891],
    [11.312785, 12.203479, 14.361582, 17.071826, 12.94894, 35.288888, 26.24005, 56.892871, 0.0, 13.007432],
    [5.89625, 11.005102, 17.345279, 22.282908, 6.36347, 28.741118, 29.904652, 51.25891, 13.007432, 0.0]
])


labels = np.array(['P11', 'P12', 'P13', 'P14','P15','P16','P17','P18','P19','P20'])

#labels = np.array(['A', 'B', 'C', 'D', 'E', 'F', 'G'])

dist_df = pd.DataFrame(dist_matrix, index=labels, columns=labels)
linkage_list = []
current_labels = {label: idx for idx, label in enumerate(labels)}

def find_closest_clusters(dist_df):
    min_dist = np.inf
    closest_clusters = None
    for i in range(len(dist_df)):
        for j in range(i + 1, len(dist_df)):
            if dist_df.iloc[i, j] < min_dist:
                min_dist = dist_df.iloc[i, j]
                closest_clusters = (dist_df.index[i], dist_df.columns[j])
    return closest_clusters, min_dist

def update_distances(dist_df, cluster1, cluster2, linkage_list, min_dist, current_labels):
    new_cluster = cluster1 + cluster2

    index1 = current_labels[cluster1]
    index2 = current_labels[cluster2]
    linkage_list.append([index1, index2, min_dist, len(new_cluster)])

    new_index = max(current_labels.values()) + 1
    current_labels[new_cluster] = new_index
    del current_labels[cluster1]
    del current_labels[cluster2]

    new_distances = []

    for cluster in dist_df.index:
        if cluster != cluster1 and cluster != cluster2:
            new_distance = min(dist_df.loc[cluster1, cluster], dist_df.loc[cluster2, cluster])
            new_distances.append(new_distance)

    new_row = pd.Series(new_distances, index=[c for c in dist_df.index if c not in [cluster1, cluster2]])
    new_row[new_cluster] = 0

    dist_df = dist_df.drop([cluster1, cluster2], axis=0).drop([cluster1, cluster2], axis=1)

    dist_df[new_cluster] = new_row
    dist_df.loc[new_cluster] = new_row.T

    return dist_df

# Realizar el clustering
with open('distancia_minima_3.txt', 'w') as f:
  while len(dist_df) > 1:
      (cluster1, cluster2), min_dist = find_closest_clusters(dist_df)
      dist_df = update_distances(dist_df, cluster1, cluster2, linkage_list, min_dist, current_labels)
      f.write(f"Clusteres: {cluster1} y {cluster2}\n")
      f.write(dist_df.to_string())
      f.write('\n\n')

linkage_matrix = np.array(linkage_list)

# Dibujar el dendrograma
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=labels, above_threshold_color='gray')

for i, d, c in zip(linkage_matrix[:, 0:2], linkage_matrix[:, 2], linkage_matrix[:, 3]):
    x = np.mean(i)
    y = d
    plt.plot(x, y, 'ro')
    plt.text(x, y, f'{y:.2f}', va='bottom', ha='center')

plt.xlabel("Clusters")
plt.ylabel("Distancia")
plt.show()

#MAXIMA

dist_df = pd.DataFrame(dist_matrix, index=labels, columns=labels)
linkage_list = []
current_labels = {label: idx for idx, label in enumerate(labels)}

def find_closest_clusters(dist_df):
    min_dist = np.inf
    closest_clusters = None
    for i in range(len(dist_df)):
        for j in range(i + 1, len(dist_df)):
            if dist_df.iloc[i, j] < min_dist:
                min_dist = dist_df.iloc[i, j]
                closest_clusters = (dist_df.index[i], dist_df.columns[j])
    return closest_clusters, min_dist

def update_distances(dist_df, cluster1, cluster2, linkage_list, min_dist, current_labels):
    new_cluster = cluster1 + cluster2

    index1 = current_labels[cluster1]
    index2 = current_labels[cluster2]
    linkage_list.append([index1, index2, min_dist, len(new_cluster)])

    new_index = max(current_labels.values()) + 1
    current_labels[new_cluster] = new_index
    del current_labels[cluster1]
    del current_labels[cluster2]

    new_distances = []

    for cluster in dist_df.index:
        if cluster != cluster1 and cluster != cluster2:
            new_distance = max(dist_df.loc[cluster1, cluster], dist_df.loc[cluster2, cluster])
            new_distances.append(new_distance)

    new_row = pd.Series(new_distances, index=[c for c in dist_df.index if c not in [cluster1, cluster2]])
    new_row[new_cluster] = 0

    dist_df = dist_df.drop([cluster1, cluster2], axis=0).drop([cluster1, cluster2], axis=1)

    dist_df[new_cluster] = new_row
    dist_df.loc[new_cluster] = new_row.T

    return dist_df

# Realizar el clustering
with open('distancia_maxima_3.txt', 'w') as f:
  while len(dist_df) > 1:
      (cluster1, cluster2), min_dist = find_closest_clusters(dist_df)
      dist_df = update_distances(dist_df, cluster1, cluster2, linkage_list, min_dist, current_labels)
      f.write(f"Clusteres: {cluster1} y {cluster2}\n")
      f.write(dist_df.to_string())
      f.write('\n\n')

linkage_matrix = np.array(linkage_list)

# Dibujar el dendrograma
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=labels, above_threshold_color='gray')

for i, d, c in zip(linkage_matrix[:, 0:2], linkage_matrix[:, 2], linkage_matrix[:, 3]):
    x = np.mean(i)
    y = d
    plt.plot(x, y, 'ro')
    plt.text(x, y, f'{y:.2f}', va='bottom', ha='center')

plt.xlabel("Clusters")
plt.ylabel("Distancia")
plt.show()

#PROMEDIO

dist_df = pd.DataFrame(dist_matrix, index=labels, columns=labels)
linkage_list = []
current_labels = {label: idx for idx, label in enumerate(labels)}

def find_closest_clusters(dist_df):
    min_dist = np.inf
    closest_clusters = None
    for i in range(len(dist_df)):
        for j in range(i + 1, len(dist_df)):
            if dist_df.iloc[i, j] < min_dist:
                min_dist = dist_df.iloc[i, j]
                closest_clusters = (dist_df.index[i], dist_df.columns[j])
    return closest_clusters, min_dist

def update_distances(dist_df, cluster1, cluster2, linkage_list, min_dist, current_labels):
    new_cluster = cluster1 + cluster2

    index1 = current_labels[cluster1]
    index2 = current_labels[cluster2]
    linkage_list.append([index1, index2, min_dist, len(new_cluster)])

    new_index = max(current_labels.values()) + 1
    current_labels[new_cluster] = new_index
    del current_labels[cluster1]
    del current_labels[cluster2]

    new_distances = []

    for cluster in dist_df.index:
        if cluster != cluster1 and cluster != cluster2:
            new_distance = (dist_df.loc[cluster1, cluster] + dist_df.loc[cluster2, cluster]) / 2
            new_distances.append(new_distance)

    new_row = pd.Series(new_distances, index=[c for c in dist_df.index if c not in [cluster1, cluster2]])
    new_row[new_cluster] = 0

    dist_df = dist_df.drop([cluster1, cluster2], axis=0).drop([cluster1, cluster2], axis=1)

    dist_df[new_cluster] = new_row
    dist_df.loc[new_cluster] = new_row.T

    return dist_df

# Realizar el clustering
with open('distancia_promedio_3.txt', 'w') as f:
  while len(dist_df) > 1:
      (cluster1, cluster2), min_dist = find_closest_clusters(dist_df)
      dist_df = update_distances(dist_df, cluster1, cluster2, linkage_list, min_dist, current_labels)
      f.write(f"Clusteres: {cluster1} y {cluster2}\n")
      f.write(dist_df.to_string())
      f.write('\n\n')

linkage_matrix = np.array(linkage_list)

# Dibujar el dendrograma
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=labels, above_threshold_color='gray')

for i, d, c in zip(linkage_matrix[:, 0:2], linkage_matrix[:, 2], linkage_matrix[:, 3]):
    x = np.mean(i)
    y = d
    plt.plot(x, y, 'ro')
    plt.text(x, y, f'{y:.2f}', va='bottom', ha='center')

plt.xlabel("Clusters")
plt.ylabel("Distancia")
plt.show()