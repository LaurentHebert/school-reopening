
% simple script that load contact matrices, modifies them, and outputs R0s
% no argument, everything is hardcoded
% @author: Laurent Hébert-Dufresne <lhebertd@uvm.edu>

% Epidemiological parameters
sigma = [0.34 0.34 0.34 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.44];
gamma = 1/5.1;

% Load contact matrices
M_Baseline = load('../data/SH_Baseline.dat' );
M_Lockdown = load('../data/SH_Lockdown.dat' );

% Set up exposure matrices
K_Baseline = 0*M_Baseline;
K_Lockdown = 0*M_Baseline;
K_LockdownWithSchools = 0*M_Baseline;
K_LockdownWithXSchools = 0*M_Baseline;
K_LockdownWithYSchools = 0*M_Baseline;
K_NoSchools = 0*M_Baseline;

% Parameters for partial school reopening
X1 = 1.0; X2 = 1/3; %Scenario X, activity level <10 and 10-19 resp.
Y = 1.0; %Scenario Y, activity level for kids under 10.

% Calculate exposure matrices
for i=1:size(M_Baseline,1)
    for j=1:size(M_Baseline,2)
        K_Baseline(i,j) = sigma(i)*M_Baseline(i,j);
        K_Lockdown(i,j) = sigma(i)*M_Lockdown(i,j);
        if(i<3 && j<3)
            K_LockdownWithSchools(i,j) = sigma(i)*M_Baseline(i,j);
            K_LockdownWithXSchools(i,j) = sigma(i)*X1*M_Baseline(i,j);
            K_LockdownWithYSchools(i,j) = sigma(i)*Y*M_Baseline(i,j);
            K_NoSchools(i,j) = sigma(i)*M_Lockdown(i,j);
        elseif(i<5 && j<5)
            K_LockdownWithSchools(i,j) = sigma(i)*M_Baseline(i,j);
            K_LockdownWithXSchools(i,j) = sigma(i)*X2*M_Baseline(i,j);
            K_LockdownWithYSchools(i,j) = sigma(i)*M_Lockdown(i,j);
            K_NoSchools(i,j) = sigma(i)*M_Lockdown(i,j);
        else
            K_LockdownWithXSchools(i,j) = sigma(i)*M_Lockdown(i,j);
            K_LockdownWithYSchools(i,j) = sigma(i)*M_Lockdown(i,j);
            K_LockdownWithSchools(i,j) = sigma(i)*M_Lockdown(i,j);
            K_NoSchools(i,j) = sigma(i)*M_Baseline(i,j);
        end
    end
end

% Spectral radius of the exposure matrices
lambda_Baseline = eigs(K_Baseline, 1, 'lm');
lambda_LockdownWithSchools = eigs(K_LockdownWithSchools, 1, 'lm')
lambda_LockdownWithXSchools = eigs(K_LockdownWithXSchools, 1, 'lm')
lambda_LockdownWithYSchools = eigs(K_LockdownWithYSchools, 1, 'lm')
lambda_Lockdown = eigs(K_Lockdown, 1, 'lm')
lambda_NoSchools = eigs(K_NoSchools, 1, 'lm')

% Calculate output
data = zeros(0,0);
for baseR0 = 0.0:0.001:6

    beta = baseR0*gamma/lambda_Baseline;

    R0_LockdownWithSchools = (beta/gamma)*lambda_LockdownWithSchools;
    R0_LockdownWithXSchools = (beta/gamma)*lambda_LockdownWithXSchools;
    R0_LockdownWithYSchools = (beta/gamma)*lambda_LockdownWithYSchools;
    R0_Lockdown = (beta/gamma)*lambda_Lockdown;
    R0_NoSchools = (beta/gamma)*lambda_NoSchools;

    data = [data; baseR0, R0_NoSchools, R0_Lockdown, R0_LockdownWithXSchools, R0_LockdownWithYSchools, R0_LockdownWithSchools];

end

% Save output
save('../data/SH_R0s.dat', 'data')

