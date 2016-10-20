function createSamples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20/10/2016
% Copyright 2016 Asad Mahmmod
% use the associated readme file to learn how to use all the matlab files
% The code is not optimized/cleaned  - please
% contact asadmehmud//at//gmail//dot//com for any related queries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc


%Importing the files required for REAL DATA
                       
                % AVIRIS JPL locator    
                     signal_real_import = importdata('.\AvCup11_smallcorrected_spec2.mat');
                    %signal_real_import = importdata('C:\asad\asad-EIE\research\simulations\hsi\data\real datasets\aviris\avcup12rad\extract2_spec4\Image_matlab_correctedwithfixed_norm_extract2_spec4.mat');
           
                     %for 2D image
                         %[npixel,bands] =size(signal_real_import);
                         %signal_real2D = signal_real_import;
                     %for 3D image
                         [rows,cols,bands] =size(signal_real_import)
                         npixel =rows*cols;
                         signal_real2D = reshape(signal_real_import,npixel,bands);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample Creation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Input Data 
    
     data_in = signal_real2D;
     data_in_cent = detrend(data_in,'constant');
    
    p=bands;
    bandnum_max_foritr=p;
    for bandnum=1:bandnum_max_foritr%p%length(bands_list_noiseest)
           
          bandnum
              
          %% initial vectors/matrices definition
                        
              data_in_del = data_in(:,setdiff([1:p],bandnum));
              data_in_del_cent = data_in_cent(:,setdiff([1:p],bandnum));
              data_in_cent_bandnum = data_in_cent(:,bandnum);
         
          

 %% Calculation of the Residuals               
 
              controlvecs_data_noconst_cent = [data_in_del_cent];
               
              %Calculating the Alphas a
                    alpha_noconst_cent_calc = controlvecs_data_noconst_cent\data_in_cent(:,bandnum);
                    alpha_noconst_cent = alpha_noconst_cent_calc;
              
              % Creating Matrices for Alphas 
                  alphamat_noconst_cent(:,bandnum) = alpha_noconst_cent;
              
             
             %Calculating the residuals (required FOR estimation of DIAGONAL TERMS in other file)
                            residual_noconst_cent_calc = data_in_cent_bandnum - (controlvecs_data_noconst_cent*alpha_noconst_cent);
                            residualmat_noconst_cent_calc(:,bandnum) = residual_noconst_cent_calc;
                
                            residualmat_noconst_cent = residualmat_noconst_cent_calc;
                            residual_noconst_cent = residualmat_noconst_cent_calc(:,bandnum);
                            
            %Calculating the residuals (required FOR estimation of OFF-DIAGONAL TERMS)
                      if(bandnum<p)
                                    bandnum_nondiag = bandnum+1;
                                    
                                    controlvecs_data_noconst_cent_nondiag = data_in_cent(:,setdiff([1:p],[bandnum, bandnum_nondiag]));

                                  %calculating the Alphas
                                        alpha_noconst_cent_calc_nondiag1 = controlvecs_data_noconst_cent_nondiag\data_in_cent(:,bandnum);
                                    
                                        alpha_noconst_cent_calc_nondiag2 = controlvecs_data_noconst_cent_nondiag\data_in_cent(:,bandnum_nondiag);
                                 
                                  
                                  alpha_noconst_cent_nondiag1 = alpha_noconst_cent_calc_nondiag1;
                                  
                                  alpha_noconst_cent_nondiag2 = alpha_noconst_cent_calc_nondiag2;
                                    

                                  %Residuals
                                        %Calculating the residuals
                                            residual_noconst_cent_calc_nondiag1 = data_in_cent(:,bandnum) - (controlvecs_data_noconst_cent_nondiag*alpha_noconst_cent_nondiag1);
                                            residual_noconst_cent_calc_nondiag2 = data_in_cent(:,bandnum_nondiag) - (controlvecs_data_noconst_cent_nondiag*alpha_noconst_cent_nondiag2);

                                            residual_noconst_cent_nondiag1 = residual_noconst_cent_calc_nondiag1;
                                            residual_noconst_cent_nondiag2 = residual_noconst_cent_calc_nondiag2;
                                            
                                  %Creating Matrices for Alphas and Residuals
                                       alphamat_noconst_cent_nondiag(:,((bandnum-1)*2)+1) = alpha_noconst_cent_nondiag1;
                                       alphamat_noconst_cent_nondiag(:,((bandnum-1)*2)+2) = alpha_noconst_cent_nondiag2;
                                     
                      
                                      residualmat_noconst_cent_nondiag(:,((bandnum-1)*2)+1) = residual_noconst_cent_nondiag1;
                                      residualmat_noconst_cent_nondiag(:,((bandnum-1)*2)+2) = residual_noconst_cent_nondiag2;
                                      
         
                      end % end of the loop to calculate the residuals for the off-diagonal term corresponding to a given diagonal term
                            
                      % norm of the residual for the diagonal elements
                          residual_noconst_cent_norm_vec(bandnum) = (1/(npixel))*((norm(residualmat_noconst_cent(:,bandnum)))^2);
                      % product of the residuals for the non-diagonal elements
                          if(bandnum<p)
                              residual_noconst_cent_prodnondiag_vec(bandnum) = (1/(npixel))*(residual_noconst_cent_nondiag1'*residual_noconst_cent_nondiag2);
                          end
            
        
               

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       
    end % end of the loop for calculation of residuals
    

    
                        
               
                        %% Saving samples
                 
                               
                                 %alpha and resid for diag elements
                                    save '.\alphamat_avcup11rad.mat' alphamat_noconst_cent;
                                    save '.\residualmat_avcup11rad.mat' residualmat_noconst_cent;
    
                                        
                                    %alpha,beta and residual mat for
                                 %non-diagonal elements (only needed for both bands removed version of our noise correlation estimation algorithm )
                                      save '.\alphamat_nondiag_avcup11rad.mat' alphamat_noconst_cent_nondiag;
                                      save '.\residualmat_nondiag_avcup11rad.mat' residualmat_noconst_cent_nondiag;
    
 

 %%      
 %%%%%%%  the end

  