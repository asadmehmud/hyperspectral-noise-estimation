
function  statsEst;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20/10/2016
% Copyright 2016 Asad Mahmood
% function to estimate the noise statistics in hypersepctral images -
% correlated noise in real datasets
% please read the associated readme to learn how to use these files
% The code is not optimized/cleaned  - please
% contact asadmehmud//at//gmail//dot//com for any related queries


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;
    % no. of diagonals of the noiose covariance matrix to be estimated (generally 2, 5 or 10)
        num_diagonals = 5;
        image_name = 'avcup11smallspec2';
 
 %% Importing the files required for this file          
          
                        
                        alphamat_noconst_cent_import = importdata('.\alphamat_avcup11rad.mat'); 
                        residualmat_noconst_cent_import = importdata('.\residualmat_avcup11rad.mat'); 
% %                 
   %% Defining parameters 
              
          
            [npixel, bands] = size(residualmat_noconst_cent_import);
            p=bands;
            bandnum_max_foritr=p;
            num_unknowns_vec(1,p)=0;
           
            
            statsourest_covmat(bands,bands) =0;
            statsourest_corrmat(bands,bands) =0;
            statsourest_mat_upper(bands,bands) =0;
            
            residest_covmat(bands,bands) =0;
            residest_corrmat(bands,bands) =0;
            residest_mat_upper(bands,bands) =0;
   
           
            num_neighbors=  num_diagonals -1
            num_unknowns=0;
                        for n_temp=1:num_diagonals
                           num_diagonal_elements_vec(n_temp) = (p-(n_temp-1));
                           num_unknowns = num_unknowns + num_diagonal_elements_vec(n_temp);
                           num_unknowns_vec(n_temp+1) = num_unknowns;
                        end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noise Stats/Correlation Estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    
  
    %% Loop for noise statistics estimation and the number of unknowns in the noise ocvariance
  
        
        %%%% preallocating some matrices for faster run
                coeffmat_stack1(num_unknowns,num_unknowns)=0;
                
                errorterm_vec_single_new(num_unknowns,1) = 0;
                        
                                                                                                                                                                                      
                residual_prodvec_single_new(num_unknowns,1) = 0;

                noise_samplestats_prodvec_new(num_unknowns,1) = 0;
                
                
        
    disp('loop for coeffmat calc')   
    for bandnum=1:bandnum_max_foritr
           
          bandnum
              
                

    
                    alpha_noconst_cent = alphamat_noconst_cent_import(:,bandnum); 
                    alpha_noconst_cent_withzero =  [alpha_noconst_cent(1:bandnum-1)' 0 alpha_noconst_cent(bandnum:((p-1)))'];
                             
                    alpha_noconst_cent_nondiag1_single_new = alpha_noconst_cent;
                    alpha_noconst_cent_nondiag1_single_new_withzero = [alpha_noconst_cent_nondiag1_single_new(1:bandnum-1)' 0  alpha_noconst_cent_nondiag1_single_new(bandnum: p-1)']; 
          
              
                     %Importing the Residuals
                             residual_noconst_cent = residualmat_noconst_cent_import(:,bandnum);
                             residual_noconst_cent_nondiag1_single_new = residual_noconst_cent;
                             residualmat_noconst_cent(:,bandnum) = residual_noconst_cent;
             
            
                            
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
             %%% Adding section for making the code generic to add any number of diagonals             
                    for bandnum_offdiag = bandnum:bandnum+(num_diagonals-1)
                                       
                                      
                                       if(bandnum_offdiag>p)
                                           continue;
                                       end
                                        % diag_num will be 1 for the
                                        % principal diagonal, 2 for the
                                        % next and so on
                                        diag_num = (bandnum_offdiag - bandnum)+1;
                                        
                                        % getting the right alpha,gamma,resid vec
                                        
                                                 alpha_noconst_cent_nondiag2_single_new = alphamat_noconst_cent_import(:,bandnum_offdiag);
                                                 alpha_noconst_cent_nondiag2_single_new_withzero = [alpha_noconst_cent_nondiag2_single_new(1:bandnum_offdiag-1)' 0  alpha_noconst_cent_nondiag2_single_new(bandnum_offdiag: p-1)']; 

                                                 residual_noconst_cent_nondiag2_single_new = residualmat_noconst_cent_import(:,bandnum_offdiag);

                                        
                                        % defining i and j for the
                                        % concerned eq. for which all the
                                        % coefficients will be calculated
                                        i = bandnum;
                                        j = bandnum_offdiag;
                                                  
                                        %%% Initializing the coeffmat with zeros
                                                   coeffmat= zeros(p,p); 
                                                   
                                                   coeffvec =[];
                                                   
                                        %%%% looping over all the concerned
                                        %%%% coefficents
                                                   for k = 1:p
                                                       %for l=k-num_neighbors:k+num_neighbors
                                                       for l=k:k+num_neighbors
                                                           k=k;
                                                           l=l;
                                                        %%%% Implementing the equation with kron function to find the 'kl'th coeff and placing it over
                                                         %%%% 'kl'th location in the coefficient matrix 
                                                           
                                                        
                                                     
                                                           if((l<1)||(l>(p)))
                                                               continue;
                                                           end
                                                           
%                                                          
                                                             % for SINGLE
                                                             % BAND removed
                                                             % case
                                                                      if(l==k)
                                                                           coeffmat(k,l)= (1)*asad_kron(i,k)*asad_kron(j,l) - ...
                                                                                                                    (alpha_noconst_cent_nondiag1_single_new_withzero(k))*asad_kron(j,l) - ...
                                                                                                                    (alpha_noconst_cent_nondiag2_single_new_withzero(k))*asad_kron(i,l) + ...
                                                                                                                    (alpha_noconst_cent_nondiag1_single_new_withzero(k))*(alpha_noconst_cent_nondiag2_single_new_withzero(l));
                                                                      end
                                                                      % for off
                                                                      % diag coeffs
                                                                        if(l>k)
                                                                        coeffmat(k,l)= (1)*asad_kron(i,k)*asad_kron(j,l) - ...
                                                                                                     (alpha_noconst_cent_nondiag1_single_new_withzero(k))*asad_kron(j,l) - ...
                                                                                                     (alpha_noconst_cent_nondiag1_single_new_withzero(l))*asad_kron(j,k) - ...
                                                                                                     (alpha_noconst_cent_nondiag2_single_new_withzero(k))*asad_kron(i,l) - ...
                                                                                                     (alpha_noconst_cent_nondiag2_single_new_withzero(l))*asad_kron(i,k) + ...
                                                                                                     (alpha_noconst_cent_nondiag1_single_new_withzero(k))*(alpha_noconst_cent_nondiag2_single_new_withzero(l)) + ...
                                                                                                     (alpha_noconst_cent_nondiag1_single_new_withzero(l))*(alpha_noconst_cent_nondiag2_single_new_withzero(k));
                                                                        end   
                                                                
                                                       end % end of for loop for l
                                                   end     % end of for loop for k  
                                                   
                                      

                                                                                                    
                                        %%% picking the principal diagonal and the adjacent one
                                        %%% and placing them in the form of a vector
                                                  for diagnum_temp = 0:num_neighbors
                                                         coeffvec_temp = diag(coeffmat,diagnum_temp);
                                                         coeffvec= [coeffvec ; coeffvec_temp];
                                                     end
                            
                                                      
                            % index refers to the row/eq number of the
                            % current coeffcient vector in the big matrix 
                                index = num_unknowns_vec(diag_num)+bandnum;                 
                             
                       %%%% Making other vectors like the residnorm vec etc. of generic length
                                                                                                                                                                                      
                                            residual_prodvec_single_new(index) = (1/(npixel))*(residual_noconst_cent_nondiag1_single_new'*residual_noconst_cent_nondiag2_single_new);

                                            %noise_samplestats_prodvec_new(index) = (1/(n*n))*(noise(:,bandnum)'*noise(:,bandnum_offdiag));

                                                      
                                                      
                                                      
                      %%%%%%%% creating matrices for the coefficient term
                      %%%%%%%% vectors %%%%%%%%

                           
                            % for orig michael's coefff eq single band
                                coeffmat_stack1(index,:) = coeffvec';
                             

                    end % end of the for loop for bandnum_offdiag       
             
             %%% end of the section for making the code generic to add any
             %%% number of diagonals
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
       
    end % end of the loop for calculation of residuals
    
    disp('calc of estimate via matrix inverse')  
    

       
  %% Estimates of the Noise Statistics (Correlated) - vector estimate via the matrix inversion
     
            
                       noise_stats_approxest = (coeffmat_stack1)\(residual_prodvec_single_new);
                    
        
%% Making a matrix from the vector estimates 
                 %create the upper traingular part of our matrix
                 statsourest_vec = noise_stats_approxest;
                 for num_diag_temp = 1: num_diagonals
                     statsourest_mat_temp_upper(bands,bands)=0;
                     residest_mat_temp_upper(bands,bands)=0;
                     
                     concerned_indices = (num_unknowns_vec(num_diag_temp)+1):num_unknowns_vec(num_diag_temp+1);
                     size_concerned_indices = size(concerned_indices);
                     concerned_diag_ourest = statsourest_vec(concerned_indices);
                     concerned_diag_residest = residual_prodvec_single_new(concerned_indices);
                     
                     statsourest_mat_temp_upper = diag(concerned_diag_ourest,num_diag_temp-1);
                     residest_mat_temp_upper = diag(concerned_diag_residest,num_diag_temp-1);
                     
                     
                     statsourest_mat_upper = statsourest_mat_upper + statsourest_mat_temp_upper; 
                     residest_mat_upper = residest_mat_upper + residest_mat_temp_upper; 
                 end
                 
                 % create a symmetric matrix
                     statsourest_mat_upper = statsourest_mat_upper;
                     statsourest_mat_lower = tril(statsourest_mat_upper',-1);
                     residest_mat_lower = tril(residest_mat_upper',-1);
                 
                     statsourest_covmat = statsourest_mat_upper+statsourest_mat_lower;
                     residest_covmat = residest_mat_upper + residest_mat_lower;
                 
                 %making correlation matrices from the covariance matrices
                     statsourest_corrmat = covtocorr(statsourest_covmat);
                     residest_corrmat = covtocorr(residest_covmat);
                     
                 % creating an image of the estimated noise correlation and
                 % covariance matrices
                 
                  set(0,'DefaultFigureColormap',feval('gray'));
                    imagesc(statsourest_covmat)
                    colorbar;
                    figure
                    set(0,'DefaultFigureColormap',feval('gray'));
                    imagesc(statsourest_corrmat)
                    colorbar;


  