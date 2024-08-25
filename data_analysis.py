import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import csv 
from scipy.spatial.distance import euclidean

# id	 x	 y	 theta	 volume	 gfp	 rfp	 yfp	 cfp
bacterias_time = []

with open('conteos.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile)
    first_row = next(spamreader)  # Read the first row (header)
    bacteria_actual_group = []
    for row in spamreader:
        if(row[0] == "id"):
            bacterias_time.append(bacteria_actual_group)
            bacteria_actual_group = []
        else:
            bacteria_actual_group.append(row)
    print("Number of times periods", len(bacterias_time))
    #print(bacterias_time[0])

# LOOK IN TIME 1 FOR BACTERIAS GFP
bacteria_gfp = []
bacteria_rfp = []
bacteria_yfp = []

for bacterium in bacterias_time[1500]:
    if(int(bacterium[5]) > 0):
        print("IS GFP!", bacterium)
        bacteria_gfp.append(bacterium)

print(bacteria_gfp)

positions = np.array([[float(bacterium[1]), float(bacterium[2])] for bacterium in bacteria_gfp])

# DBSCAN clustering
db = DBSCAN(eps=50, min_samples=2).fit(positions)

# Get cluster labels
labels = db.labels_

# Number of clusters found (-1 indicates noise points)
n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
print(f'Number of clusters: {n_clusters}')

# Organize bacteria into clusters
clusters = {}
for label, bacterium in zip(labels, bacteria_gfp):
    if label not in clusters:
        clusters[label] = []
    clusters[label].append([float(bacterium[1]), float(bacterium[2])])

# Calculate the centroid of each cluster
centroids = {}
for cluster_id, cluster_bacteria in clusters.items():
    cluster_positions = np.array(cluster_bacteria)
    centroid = cluster_positions.mean(axis=0)  # Calculate the centroid
    centroids[cluster_id] = centroid

# Print the centroids
print("\nCentroids of Clusters:")
for cluster_id, centroid in centroids.items():
    print(f"Cluster {cluster_id} Centroid: {centroid}")

# Calculate the distance between each pair of cluster centroids
distances = {}
cluster_ids = list(centroids.keys())

for i in range(len(cluster_ids)):
    for j in range(i + 1, len(cluster_ids)):
        cluster_id_1 = cluster_ids[i]
        cluster_id_2 = cluster_ids[j]
        distance = euclidean(centroids[cluster_id_1], centroids[cluster_id_2])
        distances[(cluster_id_1, cluster_id_2)] = distance

# Print the distances between clusters
print("\nDistances Between Clusters:")
for cluster_pair, distance in distances.items():
    print(f"Distance between Cluster {cluster_pair[0]} and Cluster {cluster_pair[1]}: {distance}")

# Plot the clusters with centroids and Y-axis inverted
plt.scatter(positions[:, 0], positions[:, 1], c=labels, cmap='viridis', label='Bacteria')
for cluster_id, centroid in centroids.items():
    plt.scatter(centroid[0], centroid[1], s=200, c='red', marker='X', label=f'Centroid {cluster_id}')
    plt.text(centroid[0], centroid[1], f'Centroid {cluster_id}', fontsize=12, color='red')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('DBSCAN Clustering of GFP-Positive Bacteria with Centroids')
plt.legend()
plt.gca().invert_yaxis()
plt.show()
