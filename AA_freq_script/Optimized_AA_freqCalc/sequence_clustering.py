from sklearn.cluster import KMeans
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.cluster import MeanShift
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from minisom import MiniSom
from sklearn.metrics import pairwise_distances



# k-means clustering
def kmeans_clustering(sequences, k):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=k)
    clusters = kmeans.fit_predict(amino_acid_frequencies)

    # Return the cluster assignments
    return clusters

# # Example usage:
# sequences = [
#     "ACGT",
#     "AGCT",
#     "ATCG",
#     "ACCG"
# ]
# k = 2
# clusters = kmeans_clustering(sequences, k)
# print(clusters)

# Hierarchical clustering
def hierarchical_clustering(sequences):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform hierarchical clustering
    linkage_matrix = linkage(amino_acid_frequencies, method='ward')

    # Plot the dendrogram
    plt.figure(figsize=(10, 6))
    dendrogram(linkage_matrix)
    plt.xlabel('Sequences')
    plt.ylabel('Distance')
    plt.title('Hierarchical Clustering Dendrogram')
    plt.show()

# DBSCAN (Density-Based Spatial Clustering of Applications with Noise):



def dbscan_clustering(sequences, eps, min_samples):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform DBSCAN clustering
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    clusters = dbscan.fit_predict(amino_acid_frequencies)

    # Plot the clustering result
    unique_clusters = set(clusters)
    num_clusters = len(unique_clusters)
    colors = plt.cm.get_cmap('tab10', num_clusters)

    plt.figure(figsize=(10, 6))
    for i, cluster in enumerate(unique_clusters):
        if cluster == -1:
            # Outliers
            color = 'black'
        else:
            color = colors(i)
        cluster_indices = np.where(clusters == cluster)[0]
        plt.scatter(cluster_indices, [0] * len(cluster_indices), color=color, label=f'Cluster {cluster}')

    plt.xlabel('Sequences')
    plt.ylabel('Cluster')
    plt.title('DBSCAN Clustering')
    plt.legend()
    plt.show()

# # Example usage:
# sequences = [
#     "ACGT",
#     "AGCT",
#     "ATCG",
#     "ACCG",
#     "GCGT",
#     "GCTG",
#     "GTGC"
# ]
# dbscan_clustering(sequences, eps=0.3, min_samples=2)

# Mean Shift Clustering:


def mean_shift_clustering(sequences, bandwidth):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform Mean Shift clustering
    mean_shift = MeanShift(bandwidth=bandwidth)
    clusters = mean_shift.fit_predict(amino_acid_frequencies)

    # Plot the clustering result
    unique_clusters = set(clusters)
    num_clusters = len(unique_clusters)
    colors = plt.cm.get_cmap('tab10', num_clusters)

    plt.figure(figsize=(10, 6))
    for i, cluster in enumerate(unique_clusters):
        cluster_indices = np.where(clusters == cluster)[0]
        plt.scatter(cluster_indices, [0] * len(cluster_indices), color=colors(i), label=f'Cluster {cluster}')

    plt.xlabel('Sequences')
    plt.ylabel('Cluster')
    plt.title('Mean Shift Clustering')
    plt.legend()
    plt.show()

# # Example usage:
# sequences = [
#     "ACGT",
#     "AGCT",
#     "ATCG",
#     "ACCG",
#     "GCGT",
#     "GCTG",
#     "GTGC"
# ]
# mean_shift_clustering(sequences, bandwidth=0.3)

# Affinity propagation

def affinity_propagation_clustering(sequences):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform Affinity Propagation clustering
    affinity_propagation = AffinityPropagation()
    exemplars = affinity_propagation.fit_predict(amino_acid_frequencies)

    # Plot the clustering result
    unique_exemplars = set(exemplars)
    num_exemplars = len(unique_exemplars)
    colors = plt.cm.get_cmap('tab10', num_exemplars)

    plt.figure(figsize=(10, 6))
    for i, exemplar in enumerate(unique_exemplars):
        exemplar_indices = np.where(exemplars == exemplar)[0]
        plt.scatter(exemplar_indices, [0] * len(exemplar_indices), color=colors(i), label=f'Exemplar {exemplar}')

    plt.xlabel('Sequences')
    plt.ylabel('Exemplar')
    plt.title('Affinity Propagation Clustering')
    plt.legend()
    plt.show()

# Example usage:
# sequences = [
#     "ACGT",
#     "AGCT",
#     "ATCG",
#     "ACCG",
#     "GCGT",
#     "GCTG",
#     "GTGC"
# ]
# affinity_propagation_clustering(sequences)

# spectral clustering


def spectral_clustering(sequences, n_clusters=2):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform dimensionality reduction using PCA
    pca = PCA(n_components=2)
    transformed_data = pca.fit_transform(amino_acid_frequencies)

    # Perform Spectral Clustering
    spectral_clustering = SpectralClustering(n_clusters=n_clusters)
    labels = spectral_clustering.fit_predict(transformed_data)

    # Plot the clustering result
    unique_labels = set(labels)
    colors = plt.cm.get_cmap('tab10', len(unique_labels))

    plt.figure(figsize=(8, 6))
    for i, label in enumerate(unique_labels):
        indices = np.where(labels == label)[0]
        plt.scatter(transformed_data[indices, 0], transformed_data[indices, 1], color=colors(i), label=f'Cluster {label}')

    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('Spectral Clustering')
    plt.legend()
    plt.show()

# Example usage:
# sequences = [
#     "ACGT",
#     "AGCT",
#     "ATCG",
#     "ACCG",
#     "GCGT",
#     "GCTG",
#     "GTGC"
# ]
# spectral_clustering(sequences, n_clusters=3)

# Gaussian Mixture models (GMM):


def gmm_clustering(sequences, n_components=2):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform dimensionality reduction using PCA
    pca = PCA(n_components=2)
    transformed_data = pca.fit_transform(amino_acid_frequencies)

    # Perform GMM clustering
    gmm = GaussianMixture(n_components=n_components)
    gmm.fit(amino_acid_frequencies)
    labels = gmm.predict(amino_acid_frequencies)

    # Plot the clustering result
    unique_labels = set(labels)
    colors = plt.cm.get_cmap('tab10', len(unique_labels))

    plt.figure(figsize=(8, 6))
    for i, label in enumerate(unique_labels):
        indices = np.where(labels == label)[0]
        plt.scatter(transformed_data[indices, 0], transformed_data[indices, 1], color=colors(i), label=f'Cluster {label}')

    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('Gaussian Mixture Models Clustering')
    plt.legend()
    plt.show()

# Example usage:
# sequences = [
#     "ACGT",
#     "AGCT",
#     "ATCG",
#     "ACCG",
#     "GCGT",
#     "GCTG",
#     "GTGC"
# ]
# gmm_clustering(sequences, n_components=3)

# som clustering:
def som_clustering(sequences, grid_shape=(10, 10), n_iterations=100):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Convert the list to a numpy array
    amino_acid_frequencies = np.array(amino_acid_frequencies)

    # Perform dimensionality reduction using PCA
    pca = PCA(n_components=2)
    transformed_data = pca.fit_transform(amino_acid_frequencies)

    # Initialize the SOM
    som = MiniSom(grid_shape[0], grid_shape[1], transformed_data.shape[1], sigma=1.0, learning_rate=0.5)
    som.random_weights_init(transformed_data)

    # Train the SOM
    som.train_batch(transformed_data, n_iterations)

    # Get the cluster labels for each sequence
    cluster_labels = np.zeros(transformed_data.shape[0])
    for i, sample in enumerate(transformed_data):
        bmu = som.winner(sample)
        cluster_labels[i] = bmu[1] * grid_shape[0] + bmu[0]

    # Plot the clustering result
    unique_labels = np.unique(cluster_labels)
    colors = plt.cm.get_cmap('tab10', len(unique_labels))

    plt.figure(figsize=(8, 6))
    for i, label in enumerate(unique_labels):
        indices = np.where(cluster_labels == label)[0]
        plt.scatter(transformed_data[indices, 0], transformed_data[indices, 1], color=colors(i),
                    label=f'Cluster {int(label)}')

    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('Self-Organizing Maps (SOM) Clustering')
    plt.legend()
    plt.show()


# Example usage:
# sequences = [
#     "ACGT",
#     "AGCT",
#     "ATCG",
#     "ACCG",
#     "GCGT",
#     "GCTG",
#     "GTGC"
# ]
# som_clustering(sequences, grid_shape=(5, 5), n_iterations=100)

# Affinity clustering:


def affinity_clustering(sequences, damping=0.5):
    # Calculate amino acid frequencies for each sequence
    amino_acid_frequencies = []
    max_length = max(len(seq) for seq in sequences)

    for sequence in sequences:
        amino_acid_counts = {}
        total_aa_count = 0
        for aa in sequence:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
            else:
                amino_acid_counts[aa] = 1
            total_aa_count += 1
        frequencies = [amino_acid_counts.get(aa, 0) / total_aa_count for aa in sorted(amino_acid_counts.keys())]
        # Pad the frequencies list with zeros to match the maximum length
        frequencies += [0] * (max_length - len(frequencies))
        amino_acid_frequencies.append(frequencies)

    # Calculate the pairwise affinity scores using cosine similarity
    affinity_matrix = 1 - pairwise_distances(amino_acid_frequencies, metric='cosine')

    # Perform Affinity Propagation clustering
    affinity_clustering = AffinityPropagation(damping=damping, affinity='precomputed')
    cluster_labels = affinity_clustering.fit_predict(affinity_matrix)

    # Get the unique cluster labels
    unique_labels = np.unique(cluster_labels)

    # Print the clusters
    for label in unique_labels:
        indices = np.where(cluster_labels == label)[0]
        cluster_sequences = [sequences[i] for i in indices]
        print(f'Cluster {label}:')
        for sequence in cluster_sequences:
            print(sequence)
        print()

# Example usage:
sequences = [
    "ACGT",
    "AGCT",
    "ATCG",
    "ACCG",
    "GCGT",
    "GCTG",
    "GTGC"
]
affinity_clustering(sequences, damping=0.5)

