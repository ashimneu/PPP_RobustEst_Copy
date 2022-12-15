The RAPS.m script takes following input and generates following output.
INPUT
y   - measurements
H   - measurement matrix
P   - Prior Covariance matrix
R   - Measurement Covariance matrix
J_l - Information Matrix Lower Bound
x_prior - prior state vector estimate

OUTPUT:   
x_post   - posterior state vector estimate
by       - measurement selection vector (binary)
augcost  - augmented cost for RAPS B&B
exitflag - see MATLAB function: intlinprog for description

augcost is the cost as described in eqn (11).

More details on formulation of 
ieqLHS % inequality constraint LHS matrix
ieqRHS % inequality constraint RHS vector 
are given in Optimization_inequality.pdf.

