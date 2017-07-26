####################################################
# Mariona Roigé Valiente, 2016                     #
# Ecological Informatics, BPRC, Lincoln University #              			  
# mariona.roige.valiente@gmail.com       	   #
####################################################

% This program is doing an unsupervised SOM
% The file must be in SOM_PAK format
% No target variable
% The library 'SOM_Toolbox' must be added to Matlab path
% to be able to use its functions. 

%%%%% Importing & preparing data %%%%%%%%%%%%%%%%%%%%%%

% Reading in the data
sD = som_read_data('datasetA.txt');

%Size of data
[n,p] = size(sD.data);

% Parameter selection
struct=som_topol_struct('Data.txt',sD); 

% Map initialisation
sM=som_lininit(sD,'msize',struct.msize,'lattice',struct.lattice);

% Map construction
sM=som_batchtrain(sM,sD,'trainlen',5000,'radius_ini',3,'radius_fin',1);

% Name of countries association
sM=som_autolabel(sM,sD);

% Quality map evaluation
[qe te]=som_quality(sM,sD);

% Graphic representation of the model
som_show(sM,'umat','all','empty','Labels','norm','d')

som_show_add('label',sM.labels,'textsize',6,'textcolor','r')

% Weights extraction
Bmus = som_bmus(sM,sD);
rNZ = sM.codebook(Bmus(306),:);

% Ranks computation
Rankedweights = Rank(rNZ); % Using function "Rank" made by Adrien Souquet

%Grid with the number of neuron
[n1,p1] = size(sM.codebook)
label = (1:n1)
% put into one column 
label = label(:)
label = num2cell(label)
figure
som_show(sM,'empty','Neuron','norm','d');
som_show_add('label',label,'textsize',12,'textcolor','b'),title('Grid of the number of neurons ');

