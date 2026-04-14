function Model = checkResonance(Model, Master, Slave, Indp, errorBar, style)

if nargin <= 5
    PS = 'CNF';
else
    PS = style.PS;
end

CS = Indp.type;

error = errorBar*min(abs(Model.spectrum.Lambda(Master)));
n_p = length(Indp.I);
ResonanceFlag = cell(n_p,1);
ResonanceI = cell(n_p,1);

% pi = 1的情况
p_i = 1;
l_p1 = size(Indp.I{p_i},2);
ResonanceFlag{p_i} = zeros(l_p1,1);
ResonanceI{p_i} = cell(l_p1,1);

switch CS
    case 'index'
        for p_i = 2:n_p
            Mu = change_Ind_Mu(Indp.I{p_i},length(Master));
            sumLambda = Mu'*Model.spectrum.Lambda(Master);
            l_pi = length(sumLambda);
            ResonanceFlag{p_i} = zeros(l_pi,1);
            ResonanceI{p_i}=cell(l_pi,1);
            for n_i = 1:l_pi
                R_ind = find(abs(sumLambda(n_i) - Model.spectrum.Lambda(Master))/p_i<error);
                switch PS
                    case 'GP'
                        ResonanceFlag{p_i}(n_i) = 1;
                        ResonanceI{p_i}{n_i} = 1:length(Master);
                    case 'CNF'
                        if sum(R_ind)>0
                            disp(['The inner resonance occur for',mat2str(Mu(:,n_i)),':',mat2str(R_ind),'-th master mode'])
                            ResonanceFlag{p_i}(n_i) = 1;
                            ResonanceI{p_i}{n_i} = R_ind;
                        end
                    case 'RNF'
                        if sum(R_ind)>0
                            disp(['The inner resonance occur for',mat2str(Mu(:,n_i)),':',mat2str(R_ind),'-th master mode'])
                            ResonanceFlag{p_i}(n_i) = 1;
                            if R_ind <= length(Master)/2
                                ResonanceI{p_i}{n_i} = [R_ind R_ind+length(Master)/2];
                            else
                                ResonanceI{p_i}{n_i} = [R_ind R_ind-length(Master)/2];
                            end
                        end
                end
            end
            
            for n_i = 1:l_pi
                R_ind = find(abs(sumLambda(n_i) - Model.spectrum.Lambda(Slave))<error);
                if sum(R_ind)>0
                    disp(['The outer resonance occur for',mat2str(Mu(:,n_i)),':',mat2str(R_ind),'-th slave mode'])
                end
            end
        end

    case 'multi_index'
        for p_i = 2:n_p
            Mu = Indp.I{p_i};
            sumLambda = Mu'*Model.spectrum.Lambda(Master);
            l_pi = length(sumLambda);
            ResonanceFlag{p_i} = zeros(l_pi,1);
            ResonanceI{p_i}=cell(l_pi,1);
            for n_i = 1:l_pi
                R_ind = find(abs(sumLambda(n_i) - Model.spectrum.Lambda(Master))/p_i<error);
                switch PS
                    case 'GP'
                        ResonanceFlag{p_i}(n_i) = 1;
                        ResonanceI{p_i}{n_i} = 1:length(Master);
                    case 'CNF'
                        if sum(R_ind)>0
                            disp(['The inner resonance occur for',mat2str(Mu(:,n_i)),':',mat2str(R_ind),'-th master mode'])
                            ResonanceFlag{p_i}(n_i) = 1;
                            ResonanceI{p_i}{n_i} = R_ind;
                        end
                    case 'RNF'
                        if sum(R_ind)>0
                            disp(['The inner resonance occur for',mat2str(Mu(:,n_i)),':',mat2str(R_ind),'-th master mode'])
                            ResonanceFlag{p_i}(n_i) = 1;
                            if R_ind <= length(Master)/2
                                ResonanceI{p_i}{n_i} = [R_ind R_ind+length(Master)/2];
                            else
                                ResonanceI{p_i}{n_i} = [R_ind R_ind-length(Master)/2];
                            end
                        end
                end
            end
            
            for n_i = 1:l_pi
                R_ind = find(abs(sumLambda(n_i) - Model.spectrum.Lambda(Slave))<error);
                if sum(R_ind)>0
                    disp(['The outer resonance occur for',mat2str(Mu(:,n_i)),':',mat2str(R_ind),'-th slave mode'])
                end
            end
        end
end
Indp.ResonanceFlag = ResonanceFlag;
Indp.ResonanceI = ResonanceI;
Model.Master = Master;
Model.Slave = Slave;
Model.Indp = Indp;
end