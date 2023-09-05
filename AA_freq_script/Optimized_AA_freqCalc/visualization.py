import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from weblogo import LogoData, LogoFormat
# from weblogo.seq import SeqCanvas
import logomaker
import seqlogo

f = '2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED_for_logo.fasta'

# seqlogo
def generate_sequence_logo(aligned_sequences):
    # Calculate the frequency matrix from the aligned sequences
    frequency_matrix = seqlogo.alignment_to_matrix(aligned_sequences, format='fraction')


    # Create a sequence logo using the frequency matrix
    logo = seqlogo.seqlogo(frequency_matrix)

    # Customize the appearance of the sequence logo
    logo.style_title('Sequence Logo')
    logo.style_color('classic')
    logo.style_show_consensus(True)
    logo.format_xticks(format='digit')

    # Display the sequence logo
    logo.plot()

#LogoBar
def generate_sequence_logo(aligned_sequences):
    # Create a Logo object and add the aligned sequences
    logo = logomaker.Logo(aligned_sequences)

    # Customize the appearance of the sequence logo
    logo.style_properties['color_scheme'] = 'chemistry'
    logo.style_properties['stack_aspect_ratio'] = 5
    logo.style_properties['show_spines'] = False
    logo.style_properties['show_baseline'] = False
    logo.ax.set_ylabel('Information Content')

    # Display the sequence logo
    plt.show()

# generate_sequence_logo(f)
# heatmap

def generate_heatmap(frequencies):
    # Extract the amino acids and their frequencies from the dictionary
    # amino_acids = list(amino_acid_frequencies[0].keys())
    # frequencies = np.array([list(freq.values()) for freq in amino_acid_frequencies])

    # Generate the heatmap
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(frequencies, cmap='YlOrRd')

    # Set the plot title and axis labels
    plt.title('Heatmap of Amino Acid Frequencies')
    plt.xlabel('Amino Acid')
    plt.ylabel('Position')

    # Set the x-axis ticks and labels
    ax.set_xticks(np.arange(len(amino_acids)))
    ax.set_xticklabels(amino_acids)

    # Set the y-axis ticks and labels
    ax.set_yticks(np.arange(len(frequencies)))
    ax.set_yticklabels(np.arange(1, len(frequencies) + 1))

    # Rotate the x-axis tick labels for better readability
    plt.xticks(rotation=45)

    # Create a colorbar for the heatmap
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Frequency')

    # Display the plot
    plt.tight_layout()
    plt.show()

# bar chart
def generate_bar_chart(frequencies, positions=None):
    # If positions are not provided, consider all positions in the alignment
    if positions is None:
        positions = range(len(frequencies))

    # Extract the amino acids and their frequencies from the dictionary
    # amino_acids = list(amino_acid_frequencies[0].keys())
    # frequencies = np.array([list(freq.values()) for freq in amino_acid_frequencies])

    # Filter the frequencies for the selected positions
    frequencies = frequencies[positions]

    # Generate the bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    width = 0.8 / len(positions)

    for i, pos in enumerate(positions):
        x = np.arange(len(amino_acids)) + i * width
        ax.bar(x, frequencies[i], width=width, label=f'Position {pos+1}')

    # Set the plot title and axis labels
    plt.title('Bar Chart of Amino Acid Frequencies')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')

    # Set the x-axis ticks and labels
    ax.set_xticks(np.arange(len(amino_acids)) + (width * (len(positions) - 1)) / 2)
    ax.set_xticklabels(amino_acids)

    # Set the legend and adjust the layout
    ax.legend(loc='upper right')
    plt.tight_layout()

    # Show the bar chart
    plt.show()

# stacked bar chart
def generate_stacked_bar_chart(amino_acid_frequencies):
    # Extract the amino acids and their frequencies from the dictionary
    amino_acids = list(amino_acid_frequencies[0].keys())
    frequencies = np.array([list(freq.values()) for freq in amino_acid_frequencies])

    # Generate the stacked bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    positions = np.arange(1, len(amino_acid_frequencies) + 1)
    bottom = np.zeros(len(amino_acid_frequencies))

    for i, amino_acid in enumerate(amino_acids):
        ax.bar(positions, frequencies[:, i], bottom=bottom, label=amino_acid)
        bottom += frequencies[:, i]

    # Set the plot title and axis labels
    plt.title('Stacked Bar Chart of Amino Acid Frequencies')
    plt.xlabel('Position')
    plt.ylabel('Frequency')

    # Set the legend and adjust the layout
    ax.legend(loc='upper right')
    plt.tight_layout()

    # Show the stacked bar chart
    plt.show()

# line plot
def generate_line_plot(amino_acid_frequencies, amino_acid):
    # Extract the frequencies of the specified amino acid from the dictionary
    frequencies = [freq.get(amino_acid, 0) for freq in amino_acid_frequencies]

    # Generate the line plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(range(1, len(frequencies) + 1), frequencies, marker='o')

    # Set the plot title and axis labels
    plt.title('Line Plot of Amino Acid Frequency')
    plt.xlabel('Position')
    plt.ylabel('Frequency')

    # Show the line plot
    plt.show()

# box plot
def generate_box_plot():
    # Create a list of lists containing amino acid frequencies for each position
    # data = [amino_acid_frequencies[position] for position in positions]
    positions = [1,2,3,4,5]
    td = [[12,1], [1, 2], [3, 3], [123, 4], [9, 5]]
    # Generate the box plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.boxplot(td)

    # Set the x-axis labels
    ax.set_xticklabels(positions)

    # Set the plot title and axis labels
    plt.title('Box Plot of Amino Acid Frequencies')
    plt.xlabel('Position')
    plt.ylabel('Frequency')

    # Show the box plot
    plt.show()

# generate_box_plot()
# radar plot

def generate_radar_plot(amino_acid_compositions, sequence_labels):
    """"
    The radar plot will display each amino acid as a separate axis, and the
    distance from the center of the plot indicates the frequency or abundance of that amino acid.
    Each sequence or group will have its own line on the plot, and the filled area
    represents the relative composition of amino acids.
    """
    # Extract the number of sequences and number of amino acids
    num_sequences = len(amino_acid_compositions)
    num_amino_acids = len(amino_acid_compositions[0])

    # Create angles for the radar plot
    angles = np.linspace(0, 2 * np.pi, num_amino_acids, endpoint=False).tolist()
    angles += angles[:1]  # Repeat the first angle to close the plot

    # Initialize the radar plot
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'polar': True})
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(list(range(1, num_amino_acids + 1)))
    ax.set_yticklabels([])

    # Plot each sequence/group on the radar plot
    for i in range(num_sequences):
        amino_acid_composition = amino_acid_compositions[i]
        values = amino_acid_composition + [amino_acid_composition[0]]  # Repeat the first value to close the plot
        ax.plot(angles, values, linewidth=2, label=sequence_labels[i])
        ax.fill(angles, values, alpha=0.25)

    # Add legend and title
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
    plt.title('Radar Plot of Amino Acid Composition')

    # Show the radar plot
    plt.show()


# import numpy as np
# import matplotlib.pyplot as plt

def generate_bubble_chart(amino_acid_frequencies):
    # Extract the number of positions and number of amino acids
    num_positions = len(amino_acid_frequencies)
    num_amino_acids = len(amino_acid_frequencies[0])

    # Prepare data for bubble chart
    positions = np.arange(1, num_positions + 1)
    frequencies = np.sum(amino_acid_frequencies, axis=1)
    dominant_amino_acids = np.argmax(amino_acid_frequencies, axis=1)
    color_map = ['red', 'green', 'blue', 'yellow']  # Customize color map as desired

    # Generate the bubble chart
    plt.scatter(positions, frequencies, s=frequencies*500, c=[color_map[i] for i in dominant_amino_acids],
                alpha=0.5, edgecolors='black')

    # Customize chart appearance
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.title('Bubble Chart of Amino Acid Frequencies')
    plt.grid(True)

    # Show the chart
    plt.show()

td = [[12, 1], [1, 2], [3, 3], [123, 4], [9, 5]]
generate_bubble_chart(td)

# import numpy as np
# import matplotlib.pyplot as plt
# import networkx as nx

def visualize_amino_acid_network(amino_acid_frequencies, threshold=0.1):
    # Convert the amino acid frequencies into a numpy array
    data = np.array(amino_acid_frequencies)

    # Calculate the co-occurrence or mutual information matrix
    # You can replace this step with your desired calculation method
    # For example, using numpy's correlation coefficient: correlation_matrix = np.corrcoef(data.T)
    # Or using a different method like mutual information estimation
    # The resulting matrix should be a square symmetric matrix
    correlation_matrix = np.corrcoef(data.T)

    # Create a network graph
    G = nx.Graph()

    # Add nodes to the graph
    num_amino_acids = len(amino_acid_frequencies[0])
    G.add_nodes_from(range(num_amino_acids))

    # Add edges to the graph based on the correlation matrix
    for i in range(num_amino_acids):
        for j in range(i + 1, num_amino_acids):
            if correlation_matrix[i, j] >= threshold:
                G.add_edge(i, j, weight=correlation_matrix[i, j])

    # Draw the network graph
    pos = nx.spring_layout(G)
    edge_weights = nx.get_edge_attributes(G, 'weight').values()
    node_colors = ['r' for _ in range(num_amino_acids)]  # Set node colors (e.g., red)
    node_sizes = [100 for _ in range(num_amino_acids)]  # Set node sizes
    edge_colors = ['gray' for _ in edge_weights]  # Set edge colors
    edge_widths = [2 * w for w in edge_weights]  # Set edge widths

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths)
    nx.draw_networkx_labels(G, pos)

    plt.title('Amino Acid Network')
    plt.axis('off')
    plt.show()

td = [[12, 1], [1, 2], [3, 3], [123, 4], [9, 5]]
visualize_amino_acid_network(td)


# principal component analysis (PCA plot)

def perform_pca(aligned_sequences):

    # Assuming you have a variable called 'amino_acid_frequencies' containing the amino acid frequencies
    # of aligned sequences in the format: [[seq1_aa_freq], [seq2_aa_freq], ...]

    # Call the function to generate the PCA plot
    #plot_pca(amino_acid_frequencies)

    # Perform PCA on the aligned sequences
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(aligned_sequences)

    # Plot the PCA results
    plt.scatter(pca_result[:, 0], pca_result[:, 1])
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA Plot')
    plt.show()


# Example usage:
aligned_sequences = [[0.2, 0.1, 0.3, 0.4], [0.5, 0.2, 0.1, 0.2], [0.3, 0.3, 0.2, 0.2]]
perform_pca(aligned_sequences)


def plot_pca(amino_acid_frequencies):
    # Convert the amino acid frequencies into a numpy array
    data = np.array(amino_acid_frequencies)

    # Perform PCA for dimensionality reduction
    pca = PCA(n_components=2)
    reduced_data = pca.fit_transform(data)

    # Plot the sequences in the reduced-dimensional space
    plt.scatter(reduced_data[:, 0], reduced_data[:, 1])
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA Plot of Aligned Sequences')
    plt.show()
