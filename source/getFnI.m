function FnI = getFnI(Fn, Psi, Indp, I, order, type, CS, varargin)
if nargin <= 5
    type = 'tensor';
end

if nargin <= 6
    CS = 'index';
end

if nargin >= 8
    FlagAssembly = 1;
else
    FlagAssembly = 0;
end
N = size(Psi{1}, 1);

FnI = zeros(N, 1);
Ind = Indp.I;

switch CS
    case 'index'
        p = length(I);
        if (order == 2)&&((~isempty(Fn{order}))||FlagAssembly)
            for k = 1:p-1
                p_1 = k;
                p_2 = p-k;
                I1 = I(1:k);
                I2 = I(k+1:p);
                j1 = find_index_VinA(I1,Ind{p_1});
                j2 = find_index_VinA(I2,Ind{p_2});
                if nargin < 8
                    switch type
                        case 'tensor'
                            FnIkln = double(ttv(Fn{order},{Psi{p_1}(:,j1), Psi{p_2}(:,j2)}, [2 3]));
                        case 'multi_index'
                            Index_Fn = Fn{order}.I;
                            Fn_vector = Fn{order}.vector;
                            FnIkln = zeros(N, 1);
                            for nm = 1:size(Index_Fn,2)
                                FnIkln = FnIkln + Fn_vector(:,nm)*Psi{p_1}(Index_Fn(1,nm),j1)...
                                    *Psi{p_2}(Index_Fn(2,nm),j2);
                            end
                    end
                else
                    mode = 'ELP';
                    Assembly = varargin{1};
                    projector2 = {Assembly.unconstrain_vector(Psi{p_1}(:,j1)), ...
                        Assembly.unconstrain_vector(Psi{p_2}(:,j2))};
                    FnIkln = Assembly.constrain_vector(Assembly.tensor_projection_within2ton('T2_project_optimized',...
                        projector2,size(projector2{1},1), [2 3], mode));
                end
                FnI = FnI + FnIkln;
            end

        elseif (order == 3)&&((~isempty(Fn{order}))||FlagAssembly)
             for k = 1:p-2
                 for l = 1:p-k-1
                    p_1 = k;
                    p_2 = l;
                    p_3 = p-k-l;
        
                    I1 = I(1:k);
                    I2 = I(k+1:k+l);
                    I3 = I(k+l+1:p);
                    j1 = find_index_VinA(I1,Ind{p_1});
                    j2 = find_index_VinA(I2,Ind{p_2});
                    j3 = find_index_VinA(I3,Ind{p_3});
                    
                    if nargin < 8
                        switch type
                            case 'tensor'
                                FnIkln = double(ttv(Fn{order},{Psi{p_1}(:,j1), Psi{p_2}(:,j2), Psi{p_3}(:,j3)}, [2 3 4]));
                            case 'multi_index'
                                Index_Fn = Fn{order}.I;
                                Fn_vector = Fn{order}.vector;
                                FnIkln = zeros(N, 1);
                                for nm = 1:size(Index_Fn,2)
                                    FnIkln = FnIkln + Fn_vector(:,nm)*Psi{p_1}(Index_Fn(1,nm),j1)...
                                        *Psi{p_2}(Index_Fn(2,nm),j2)*Psi{p_3}(Index_Fn(3,nm),j3);
                                end
                        end
                    else
                        mode = 'ELP';
                        Assembly = varargin{1};
                        projector3 = {Assembly.unconstrain_vector(Psi{p_1}(:,j1)), ...
                            Assembly.unconstrain_vector(Psi{p_2}(:,j2)), ...
                            Assembly.unconstrain_vector(Psi{p_3}(:,j3))};
                        FnIkln = Assembly.constrain_vector(Assembly.tensor_projection_within2ton('T3_project_optimized',...
                            projector3,size(projector3{1},1), [2 3 4], mode));
                    end
                    FnI = FnI + FnIkln;
                 end
            end
        end

    case 'multi_index'
        p = sum(I);
        if (order == 2)&&((~isempty(Fn{order}))||FlagAssembly)
            for k = 1:p-1
                p_1 = k;
                p_2 = p-k;

                mp1 = size(Ind{p_1},2);
                for k1 = 1:mp1
                    I1 = Ind{p_1}(:,k1);
                    I2 = I - I1;
                    if min(I2) >= 0
                    j1 = find_index_VinA(I1,Ind{p_1});
                    j2 = find_index_VinA(I2,Ind{p_2});
                    if nargin < 8
                        switch type
                            case 'tensor'
                                FnIkln = double(ttv(Fn{order},{Psi{p_1}(:,j1), Psi{p_2}(:,j2)}, [2 3]));
                            case 'multi_index'
                                Index_Fn = Fn{order}.I;
                                Fn_vector = Fn{order}.vector;
                                FnIkln = zeros(N, 1);
                                for nm = 1:size(Index_Fn,2)
                                    FnIkln = FnIkln + Fn_vector(:,nm)*Psi{p_1}(Index_Fn(1,nm),j1)...
                                        *Psi{p_2}(Index_Fn(2,nm),j2);
                                end
                        end
                    else
                        mode = 'ELP';
                        Assembly = varargin{1};
                        projector2 = {Assembly.unconstrain_vector(Psi{p_1}(:,j1)), ...
                            Assembly.unconstrain_vector(Psi{p_2}(:,j2))};
                        FnIkln = Assembly.constrain_vector(Assembly.tensor_projection_within2ton('T2_project_optimized',...
                            projector2,size(projector2{1},1), [2 3], mode));
                    end
                    FnI = FnI + FnIkln;
                    else
                        continue
                    end
                end
            end

        elseif (order == 3)&&((~isempty(Fn{order}))||FlagAssembly)
             for k = 1:p-2
                 for l = 1:p-k-1
                    p_1 = k;
                    p_2 = l;
                    p_3 = p-k-l;
                    
                    mp1 = size(Ind{p_1},2);
                    mp2 = size(Ind{p_2},2);
                    for k1 = 1:mp1
                        for k2 = 1:mp2
                            I1 = Ind{p_1}(:,k1);
                            I2 = Ind{p_2}(:,k2);
                            I3 = I - I1 - I2;
                            if min(I3) >= 0
                                j1 = find_index_VinA(I1,Ind{p_1});
                                j2 = find_index_VinA(I2,Ind{p_2});
                                j3 = find_index_VinA(I3,Ind{p_3});
                                
                                if nargin < 8
                                    switch type
                                        case 'tensor'
                                            FnIkln = double(ttv(Fn{order},{Psi{p_1}(:,j1), Psi{p_2}(:,j2), Psi{p_3}(:,j3)}, [2 3 4]));
                                        case 'multi_index'
                                            Index_Fn = Fn{order}.I;
                                            Fn_vector = Fn{order}.vector;
                                            FnIkln = zeros(N, 1);
                                            for nm = 1:size(Index_Fn,2)
                                                FnIkln = FnIkln + Fn_vector(:,nm)*Psi{p_1}(Index_Fn(1,nm),j1)...
                                                    *Psi{p_2}(Index_Fn(2,nm),j2)*Psi{p_3}(Index_Fn(3,nm),j3);
                                            end
                                    end
                                else
                                    mode = 'ELP';
                                    Assembly = varargin{1};
                                    projector3 = {Assembly.unconstrain_vector(Psi{p_1}(:,j1)), ...
                                        Assembly.unconstrain_vector(Psi{p_2}(:,j2)), ...
                                        Assembly.unconstrain_vector(Psi{p_3}(:,j3))};
                                    FnIkln = Assembly.constrain_vector(Assembly.tensor_projection_within2ton('T3_project_optimized',...
                                        projector3,size(projector3{1},1), [2 3 4], mode));
                                end
                                FnI = FnI + FnIkln;
                            else
                                continue
                            end
                        end
                    end
                 end
            end
       end
 end
end