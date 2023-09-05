from flask import Flask, render_template, request

app = Flask(__name__)

# Add your functions for calculating amino acid frequencies here


# Define the route for the home page
@app.route('/')
def home():
    return render_template('index.html')

# Define the route to handle the form submission and display results
@app.route('/calculate', methods=['POST'])
def calculate():
    # Retrieve user inputs from the form
    sequence = request.form.get('sequence', '')
    # Add other inputs as needed for your application

    # Call your function to calculate amino acid frequencies
    # Replace the example result with the actual calculated frequencies
    amino_acid_frequencies = {'A': 0.1, 'C': 0.2, 'D': 0.05, 'E': 0.15}

    return render_template('result.html', sequence=sequence, frequencies=amino_acid_frequencies)

if __name__ == '__main__':
    app.run(debug=True)

