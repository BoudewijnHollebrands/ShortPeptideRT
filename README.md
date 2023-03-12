# ShortPeptideRT
predict RT of short peptides

This is an R-code for a short peptide retention time (RT) prediction model. The code starts with some meta-data, including the title, author, date, version, and the address where the code can be found.

The model takes a set of peptide sequences and predicts the RT for each sequence using a pre-trained model. The model works as follows:

A set of non-sequence specific descriptors are calculated for each peptide. These descriptors include the length of the peptide, molecular weight, LogP (a measure of hydrophobicity), and isoelectric point (pI).

The ASP descriptors are calculated. ASP descriptors are descriptors based on the amino acid sequence of the peptide. They describe the average values of certain properties for the amino acids at the N-terminal, C-terminal, and middle positions of the peptide. The ASP descriptors calculated in this model are based on two amino acid indices obtained from the AAindexDB database.

The descriptors are used as input for a pre-trained support vector regression model (fit_svm). This model was trained on a set of peptides with known RT values.

The code then goes on to calculate the descriptors for the input peptides and predict their RT values using the pre-trained model.
