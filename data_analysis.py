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


# Extracting the spatial positions (x, y) as features for clustering
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

# List of centroid pairs to connect with lines
pairs = [
    (4, 2), (4, 3), (3, 2), (3, 6), (3, 5), (2, 1),
    (2, 5), (5, 6), (1, 5), (5, 0), (5, 1), (1, 0)  # Added (1, 0) to the list
]

# Plot the clusters with centroids and Y-axis inverted
plt.scatter(positions[:, 0], positions[:, 1], c=labels, cmap='viridis', label='Bacteria')
for cluster_id, centroid in centroids.items():
    plt.scatter(centroid[0], centroid[1], s=200, c='red', marker='X', label=f'Centroid {cluster_id}')
    plt.text(centroid[0], centroid[1], f'Centroid {cluster_id}', fontsize=12, color='red')

# Draw lines between specified pairs of centroids
for cluster_id_1, cluster_id_2 in pairs:
    centroid_1 = centroids.get(cluster_id_1)
    centroid_2 = centroids.get(cluster_id_2)
    if centroid_1 is not None and centroid_2 is not None:
        plt.plot([centroid_1[0], centroid_2[0]], [centroid_1[1], centroid_2[1]], 'k--', lw=2)  # 'k--' is a dashed black line
        distance = euclidean(centroid_1, centroid_2)
        plt.text((centroid_1[0] + centroid_2[0]) / 2, (centroid_1[1] + centroid_2[1]) / 2, 
                 f'{distance:.2f}', fontsize=12, color='blue', ha='center')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('DBSCAN Clustering of GFP-Positive Bacteria with Centroids and Lines')
plt.legend()
plt.gca().invert_yaxis()
plt.show()