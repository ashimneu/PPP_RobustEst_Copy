function p = initOutlierparam(p)
if p.eb_outlier == 1
    
    p.multisim = 1;
    % Initialize parameters for measurement outliers.
    p.outlierdbpath = './PVA/outliers/outliersdb1.mat';    
    
    % Single simulation (ordinary simulation)
    % generation mode 2
    if p.genOutlier == 1
        % generate & add outliers to measurement
        p.outlierparam.mean  = 0.5;
        p.outlierparam.width = 0.1;
        p.outlierparam.count = 3;    
    else
        % use outliers from a database file. Note: Any outlier databse should
        % only be used with its corresponding GNSS measurement database
        % it was generated for.
        if exist(p.outlierdbpath,'file')==2 % Check if the database already exists
            outlierdb = load(p.outlierdbpath);
        else
            outlierdb = [];
        end
        if ~isempty(outlierdb)
            p.outliervec   = outlierdb.vec;
            p.outlierbin   = outlierdb.bin;
            p.outlierparam = outlierdb.param;            
        else
            error('Outlier database file not found! \nPlease provide a correct path to the database.')
        end
    end
    
    % Multi simulation outlier parameter list 
    p.outlierparam.meanlist  = 0.5.*[0.5,1,2,3,4,5]; %#ok<*NBRAK>
    p.outlierparam.widthlist = 0.1.*[1:5];
    p.outlierparam.countlist = 1:8;    
    
end 
end