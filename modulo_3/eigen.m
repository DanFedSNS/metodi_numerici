clear;
Nt = 120;
tau = (Nt/2 - 1) - (Nt/4 + 2);
filename = sprintf('Nt%d_simbeta10.0.dat', Nt);

cd('C:\Users\marco\Documents\metodi_numerici\modulo_3\analysis');
fid = fopen(filename, 'r');
if fid == -1
    error('File could not be opened');
end

header = fgetl(fid);
data = fscanf(fid, '%f %f', [2, Inf]); % Read as a 2-row matrix
fclose(fid);
data = data';
first_column = data(:, 1);


x2 = first_column(2);
x4 = first_column(4);
x6 = first_column(5);

A = [first_column(6), 0, first_column(7);
     0, first_column(8) - x2^2, 0;
     first_column(7), 0, first_column(9)];


B = [first_column(6+4), 0, first_column(7+4);
     0, first_column(8+4) - x2^2, 0;
     first_column(7+4), 0, first_column(9+4)];


[eigenvectors, eigenvalues_matrix] = eig(A, B);

eigenvalues = diag(eigenvalues_matrix);

disp(A);
disp(B);
disp(eigenvalues);
gaps = -log(eigenvalues)/ tau;
fprintf("The gaps are (gap0 = %f):\n", min(gaps));
disp(gaps/min(gaps));


