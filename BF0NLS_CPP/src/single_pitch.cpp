#include "single_pitch.hpp"

single_pitch::single_pitch(int L, int N, FTYPE * pB, double trans_std_var, FTYPE voicing_prob){
    
    maxModelOrder = L;
    nFftGrid = pow(2, 14);
    //   nFftGrid= pow(2, ceil(log(nFftGrid)/log(2)));
    int F=nFftGrid;
    nData = N;
    
    pitchBounds[0] = pB[0];
    pitchBounds[1] = pB[1];
    
    int iOrder;
    int nPitchesOld=-1;
    int nPitchesOldOld=-1;
    
    // Compute some values for later
    minFftIndex = (int)round(nFftGrid*pitchBounds[0]);
    int maxFftIndex = (int)floor(nFftGrid*pitchBounds[1]);
    int minPitchGrid = minFftIndex;
    int maxPitchGrid = maxFftIndex;
    
    int nPitches = maxFftIndex - minFftIndex + 1;
    
    omega_0h = new FTYPE [maxModelOrder];
    lnBF = new FTYPE [maxModelOrder+1]; //include the NULL order
    nPitchesAll = new int [maxModelOrder];
    M = nPitches;
    Mp = nFftGrid/2+1;
    
    /* allocate */
    Gamma1 = vector(maxModelOrder*Mp); //Could be smaller
    Gamma2 = vector(maxModelOrder*Mp); //Could be smaller
    
    lsSol1 = vector(Mp);    //Could be smaller
    lsSol2 = vector(Mp);    //Could be smaller
    // // // // // // // // // // // // // // // // // // // // // // // // //
    //   added by Liming Shi
//    ArrayXXf scaled_alpha_buffer_eigen=ArrayXXf::Zero(maxModelOrder,M);
//    ArrayXXf scaled_alpha_buffer2_eigen=ArrayXXf::Ones(maxModelOrder,M);
//
//    scaled_alpha_buffer2_eigen=scaled_alpha_buffer2_eigen/maxModelOrder/M;
////    std::cout << scaled_alpha_buffer2_eigen << std::endl;
//    MatrixXd A_eigen=MatrixXd::Ones(M,M);
////     A_eigen.unaryExpr(ptr_fun(ramp))
    
    
    scaled_alpha_buffer=new FTYPE*[maxModelOrder];
    scaled_alpha_buffer2=new FTYPE*[maxModelOrder];
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        scaled_alpha_buffer[ii-1]=vector(M);
         scaled_alpha_buffer2[ii-1]=vector(M);
    }
    double prior_pitch_order=1.0/maxModelOrder/M;
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        for( int jj = 1 ; jj <= M ; ++jj){
            scaled_alpha_buffer2[ii-1][jj-1]=prior_pitch_order;
        }
    }
    scaled_alpha_buffer[0][0]=0.0/0.0;
    unvoicing_scaled_alpha_buffer=0.0/0.0;

    //     scaled_alpha_buffer[0]=0;
    logpi = vector(M);
    state_prior=vector(M);
    for( int ii = 1 ; ii <= M ; ++ii){
        logpi[ii-1]=log((double)voicing_prob/double(M)/double(maxModelOrder));
    }
    //    generating the Transition matrix for Pitch
    A = new FTYPE*[M];
    for( int ii = 1 ; ii <= M ; ++ii){
        A[ii-1]=vector(M);
    }
    temp_normalization=vector(M);
    for( int ii = 1 ; ii <= M ; ++ii){
        for( int jj = 1 ; jj <= M ; ++jj){
            A[ii-1][jj-1]=normpdf((double)(minPitchGrid+jj-1)/(double)F,(double)(minPitchGrid+ii-1)/(double)F,trans_std_var);
        }
        temp_normalization[ii-1]= sum_vec( A[ii-1],M);
    }
    for( int ii = 1 ; ii <= M ; ++ii){
        for( int jj = 1 ; jj <= M ; ++jj){
            A[ii-1][jj-1]=A[ii-1][jj-1]/temp_normalization[ii-1];
        }
    }
    
    //    generating the Transition matrix for order
    B = new FTYPE*[maxModelOrder];
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        B[ii-1]=vector(maxModelOrder);
    }
    temp_normalization=vector(maxModelOrder);
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        for( int jj = 1 ; jj <= maxModelOrder ; ++jj){
            B[ii-1][jj-1]=normpdf((double)(1+jj-1),(double)(1+ii-1),1.0);
        }
        temp_normalization[ii-1]= sum_vec( B[ii-1],maxModelOrder);
    }
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        for( int jj = 1 ; jj <= maxModelOrder ; ++jj){
            B[ii-1][jj-1]=B[ii-1][jj-1]/temp_normalization[ii-1];
        }
    }
    //    generating the Transition matrix for Voicing
    C = new FTYPE*[2];
    for( int ii = 1 ; ii <= 2 ; ++ii){
        C[ii-1]=vector(2);
    }
    C[0][0]=0.7;
    C[0][1]=1-C[0][0];
    C[1][1]=0.6;
    C[1][0]=1-C[1][1];
    
    returned_value=vector(3);

    
    logModelPrior=vector(2);
    logModelPrior[0]=voicing_prob;
    logModelPrior[1]=1-voicing_prob;

    // // // // // // // // // // // // // // // // // // // // // // // // //
    
#ifdef MKL
    t1 = vector(2*Mp); //can be used as a complex array of length Mp
    t2 = vector(Mp);
#else
    t1 = vector(Mp);   //Could be smaller
    t2 = vector(Mp);   //Could be smaller
#endif
    t3 = vector(Mp);   //Could be smaller
    t4 = vector(Mp);   //Could be smaller
    
    range = vector(Mp);
    for( int k = 0 ; k < Mp; ++k )
        range[k] = (FTYPE) k;
    
    /* Prepare the dft */
    xp = vector(F);
    CvmdInit(F, 0.0, xp); // make everything zero for zero padding later
    dftData =  (FTYPE*)FFTW_MALLOC(sizeof(FFTW_COMPLEX)*F);
    p1 = FFTW_PLAN_DFT_R2C_1D(F, xp, (FFTW_COMPLEX*)dftData, FFTW_MEASURE);
    
    /*
     fftShiftVector = ...
     exp(1i*2*pi*(0:ceil(nFftGrid/2)-1)'*(nData-1)/(2*nFftGrid));
     dftData = X(1:ceil(nFftGrid/2)).*fftShiftVector;
     */
    // just using Mp = nFftGrid/2 + 1 would be sufficient
    fftShiftVector = vector(2*Mp);
    
    CvmdExpRange(nFftGrid/2+1, range,
                 (2*M_PI*(nData-1))/(2*nFftGrid),
                 t1, t2, fftShiftVector);
    
    dftData1  = vector(Mp);
    dftData2 = vector(Mp);
    
    /* Compute crossCorrelationVectors */
    compute_cc();
    
    /* precompute gamma */
    FTYPE * phi1 = vector(Mp);    // nPitches at iOrder = 1 is sufficient
    FTYPE * psi1 = vector(Mp);    //Could be smaller
    FTYPE * alpha1 = vector(Mp);  //Could be smaller
    FTYPE * alpha2 = vector(Mp);  //Could be smaller
    FTYPE * lambda1 = vector(Mp); //Could be smaller
    FTYPE * mu1 = vector(Mp);     //Could be smaller
    
    
    costFunctionMatrix = new FTYPE*[maxModelOrder];
    //      added by Liming Shi
    bar_alpha  = new FTYPE*[maxModelOrder];
    
    for( iOrder = 1 ; iOrder <= maxModelOrder ; ++iOrder){
        //precompute some quantities
        maxFftIndexOld = maxFftIndex;
        maxFftIndex = (int)floor(MIN(nFftGrid*pitchBounds[1],
                                     ((FTYPE)nFftGrid)/(2*iOrder)-1));
        
        nPitchesOldOld = nPitchesOld;
        nPitchesOld = nPitches;
        nPitches = nPitches - (maxFftIndexOld - maxFftIndex);
        nPitchesAll[iOrder-1] = nPitches;
        
        maxPitchGrid = minPitchGrid + nPitches -1;
        
        //allocate for special sized structure
        costFunctionMatrix[iOrder-1] = vector(nPitches);
        //      added by Liming Shi
        bar_alpha[iOrder-1]=vector(M);
        
        
        
        update_gamma(iOrder, nPitches, nPitchesOld, nPitchesOldOld,
                     phi1, psi1, alpha1, lambda1,
                     mu1, Gamma1);
        
        update_gamma_p(iOrder, nPitches, nPitchesOld, nPitchesOldOld,
                       alpha2, Gamma2);
        
    }//end of iOrder for-loop
    
    del_vector(phi1);
    del_vector(psi1);
    
    del_vector(alpha1);
    del_vector(alpha2);
    del_vector(lambda1);
    del_vector(mu1);
    
}


void single_pitch::compute_cc(){
    
    int maxFftIndex = (int)floor(nFftGrid*pitchBounds[1]);
    int minPitchGrid = minFftIndex;
    int maxPitchGrid = maxFftIndex;
    
    int nPitches = maxFftIndex - minFftIndex + 1;
    
    int iOrder = 0;
    int k;
    
    crossCorrelationVectors = new FTYPE*[2*(maxModelOrder+1)+1];
    
    crossCorrelationVectors[0] = vector(nPitches);
    CvmdInit(nPitches, (1+T0_REG)*nData*0.5, crossCorrelationVectors[0]);
    
    for( int k = 1 ; k < 2*(maxModelOrder+1) + 1; k += 2 ){
        //precompute some quantities
        iOrder++;
        int maxFftIndexOld = maxFftIndex;
        maxFftIndex = (int)floor(MIN(nFftGrid*pitchBounds[1],
                                     ((FTYPE)nFftGrid)/(2*iOrder)-1));
        
        nPitches = nPitches - (maxFftIndexOld - maxFftIndex);
        
        maxPitchGrid = minPitchGrid + nPitches -1;
        
        crossCorrelationVectors[k] = vector(nPitches);
#ifdef MKL
        CvmdSinRange(nPitches, range+minPitchGrid,
                     (M_PI*nData*k)/nFftGrid, t4, t1);
        CvmdSinRange(nPitches, range+minPitchGrid,
                     (M_PI*k)/nFftGrid, t4, t2);
#else
        CvmdSinRange(minPitchGrid, maxPitchGrid,
                     (M_PI*nData*k)/nFftGrid, t1);
        CvmdSinRange(minPitchGrid, maxPitchGrid,
                     (M_PI*k)/nFftGrid, t2);
#endif
        
        
        CvmdScal(nPitches, 0.5, t1);
        CvmdDiv(nPitches, t1, t2, crossCorrelationVectors[k]);
        
        crossCorrelationVectors[k+1] = vector(nPitches);
#ifdef MKL
        CvmdSinRange(nPitches, range+minPitchGrid,
                     (M_PI*nData*(k+1))/nFftGrid, t4, t1);
        CvmdSinRange(nPitches, range+minPitchGrid,
                     (M_PI*(k+1))/nFftGrid, t4, t2);
#else
        CvmdSinRange(minPitchGrid, maxPitchGrid,
                     (M_PI*nData*(k+1))/nFftGrid, t1);
        CvmdSinRange(minPitchGrid, maxPitchGrid,
                     (M_PI*(k+1))/nFftGrid, t2);
#endif
        CvmdScal(nPitches, 0.5, t1);
        CvmdDiv(nPitches, t1, t2, crossCorrelationVectors[k+1]);
        
    }
    
    hf = new FTYPE*[2*(maxModelOrder)];
    tf = new FTYPE*[2*(maxModelOrder)];
    
    for( k = 0; k < 2*maxModelOrder ; ++k){
        tf[k] = crossCorrelationVectors[k];
        hf[k] = crossCorrelationVectors[k+2];
    }
    
}


void single_pitch::compute(FTYPE * x){
    
    /* Initialize */
    int nPitchesOld=-1;
    int iOrder;
    
    int maxFftIndex = (int)floor(nFftGrid*pitchBounds[1]);
    int minPitchGrid = minFftIndex;
    int maxPitchGrid = maxFftIndex;
    int nPitches = maxFftIndex - minFftIndex + 1;
    
    /* The dft of the data */
    CvmdCopy(nData, x, xp);
    FFTW_EXECUTE(p1); //the result is in dftData
    
    CvmdMulz(nFftGrid/2+1, fftShiftVector, dftData, t1, dftData1, dftData2);
    //dftData1 and dftData2 now contains the real and imaginary part
    //of the time shifted dft coefficients (phase rotation in frequency)
    
    
    for( iOrder = 1 ; iOrder <= maxModelOrder ; ++iOrder ){
        
        //Update the considered frequency range
        maxFftIndexOld = maxFftIndex;
        maxFftIndex = (int)floor(MIN(nFftGrid*pitchBounds[1],
                                     ((FTYPE)nFftGrid)/(2*iOrder)-1));
        
        nPitchesOld = nPitches;
        nPitches = nPitches - (maxFftIndexOld - maxFftIndex);
        
        maxPitchGrid = minPitchGrid + nPitches -1;
        
        /*
         dftDataMatrix1(iOrder, validPitchIndices) = ...
         real(dftData(validFftIndices*iOrder+1))';
         dftDataMatrix2(iOrder, validPitchIndices) = ...
         imag(dftData(validFftIndices*iOrder+1))';
         No longer saves dftDataMatrix1/dftDataMatrix2 explicit
         just get pulled out of dftData1 and dftData2
         */
        
        if( iOrder == 1 ){
            /*
             lsSol1 = dftDataMatrix1(iOrder, validPitchIndices).*gamma1;
             lsSol2 = dftDataMatrix2(iOrder, validPitchIndices).*gamma2;
             */
            CvmdPack(nPitches, dftData1 + minFftIndex*iOrder,
                     iOrder, t3);
            CvmdMul(nPitches, t3, Gamma1, lsSol1);
            
            CvmdPack(nPitches, dftData2 + minFftIndex*iOrder,
                     iOrder, t3);
            CvmdMul(nPitches, t3, Gamma2, lsSol2);
        }
        else{
            CvmdPack(nPitches, dftData1 + minFftIndex*iOrder,
                     iOrder, t3); //update_ls_sol does not use t3
            update_ls_sol(iOrder, nPitches, nPitchesOld, true,
                          t3, Gamma1, lsSol1);
            
            CvmdPack(nPitches, dftData2 + minFftIndex*iOrder,
                     iOrder, t3); //update_ls_sol does not use t3
            update_ls_sol(iOrder, nPitches, nPitchesOld, false,
                          t3, Gamma2, lsSol2);
            
        }
        
        /* Compute the cost function
         costFunctionMatrix(iOrder,validPitchIndices) = ...
         sum([dftDataMatrix1(1:iOrder, validPitchIndices);...
         dftDataMatrix2(1:iOrder, validPitchIndices)].* ...
         [lsSol1;lsSol2], 1);
         */
        
        CvmdInit(nPitches, 0.0, t4);
        for( k = 0 ; k < iOrder  ; k++){
            CvmdPack(nPitches, dftData1 + minFftIndex*(k+1), k+1, t2);
            CvmdMul(nPitches, t2, lsSol1 + k*nPitches, t1);
            CvmdAddInplace(nPitches, t4, t1);
            
            CvmdPack(nPitches, dftData2 + minFftIndex*(k+1), k+1, t2);
            CvmdMul(nPitches, t2, lsSol2 + k*nPitches, t1);
            CvmdAddInplace(nPitches, t4, t1);
        }
        CvmdCopy(nPitches, t4, costFunctionMatrix[iOrder-1]);
    }
}


single_pitch::~single_pitch(){
    
    delete [] nPitchesAll;
    
    del_vector(dftData1);
    del_vector(dftData2);
    
    del_vector(t1);
    del_vector(t2);
    del_vector(t3);
    del_vector(t4);
    
    del_vector(xp);
    FFTW_DESTROY_PLAN(p1);
    FFTW_FREE(dftData);
    
    for(int k = 0 ; k < 2*(maxModelOrder +1) + 1; ++k){
        del_vector(crossCorrelationVectors[k]);
    }
    delete [] crossCorrelationVectors;
    
    delete [] hf;
    delete [] tf;
    
    del_vector(Gamma1);
    del_vector(Gamma2);
    
    for(int k = 0 ; k < maxModelOrder; ++k){
        del_vector(costFunctionMatrix[k]);
    }
    delete [] costFunctionMatrix;
    
    del_vector(lsSol1);
    del_vector(lsSol2);
    
    del_vector(fftShiftVector);
    
    delete [] omega_0h;
    delete [] lnBF;
    
    del_vector(range);
    //   added by Liming
    del_vector(logpi);
    del_vector(temp_normalization);
    delete [] scaled_alpha_buffer;
    delete [] scaled_alpha_buffer2;
    del_vector(state_prior);
    del_vector(returned_value);
    del_vector(logModelPrior);
    for(int k = 0 ; k < M; ++k){
        del_vector(A[k]);
    }
    delete [] A;
    delete [] B;
    delete [] C;
    for(int k = 0 ; k < maxModelOrder; ++k){
        del_vector(bar_alpha[k]);
    }
    delete [] bar_alpha;
}


void single_pitch::update_ls_sol(int iOrder, int nPitches, int nPitchesOld, 
                                 bool add, FTYPE * dftData,
                                 FTYPE * Gamma, FTYPE * lsSol){
    
    /*
     lambda= (ones(iOrder,1)*(dftDataMatrix(iOrder, validPitchIndices)-...
     sum(R(1:end-1,:).*lsSol(:, validPitchIndices),1)));
     
     lsSol = [lsSol(:,validPitchIndices);zeros(1,nPitches)]...
     + lambda.*gamma;
     */
    
    CvmdInit(nPitches, 0.0, t4);
    for( k = 0 ; k < iOrder-1 ; ++k ){
        
        if(add)
            CvmdAdd(nPitches, tf[iOrder-1-k], hf[k+iOrder-1], t1);
        else
            CvmdSub(nPitches, tf[iOrder-1-k], hf[k+iOrder-1], t1);
        
        CvmdMul(nPitches, t1, lsSol + k*nPitchesOld, t2);
        CvmdAddInplace(nPitches, t4, t2);
    }
    CvmdSub(nPitches, dftData, t4, t1);
    //lambda rep is now in t1
    
    CvmdCopy((iOrder-1)*nPitchesOld, lsSol, t4); //could copy less
    for( k = 0 ; k < iOrder - 1 ; ++k ){
        CvmdMul(nPitches, t1, Gamma + Mp*(iOrder-1) + k*nPitches, t2);
        CvmdAdd(nPitches, t4 + k*nPitchesOld, t2, lsSol + k*nPitches);
    }
    CvmdMul(nPitches, t1, Gamma + Mp*(iOrder-1) + (iOrder-1)*nPitches,
            lsSol + (iOrder-1)*nPitches);
}


void single_pitch::update_gamma(int iOrder, int nPitches, int nPitchesOld, 
                                int nPitchesOldOld, FTYPE * phi,
                                FTYPE * psi, FTYPE * alpha, FTYPE * lambda,
                                FTYPE * mu, FTYPE * Gamma){
    
    if( iOrder == 1 ){
        CvmdAdd(nPitches, tf[0], hf[0], t1);
        CvmdInverse(nPitches, t1, Gamma);
        CvmdCopy(nPitches, Gamma, psi);
        CvmdMul(nPitches, tf[1], Gamma, phi);
    }
    else if( iOrder == 2 ){
        CvmdAdd(nPitches, tf[0], hf[0], t1); // R00 = t[0] + h[0];
        CvmdAdd(nPitches, tf[1], hf[1], t2); // R01 = t[1] + h[1];
        CvmdAdd(nPitches, tf[0], hf[2], t3); // R11 = t[0] + h[2];
        CvmdMul(nPitches, t1, t3, t3);
        CvmdMul(nPitches, t2, t2, t4);
        CvmdSubInplace(nPitches, t3, t4);
        CvmdInverse(nPitches, t3, t4);
        
        CvmdMul(nPitches, t2, Gamma, alpha);
        
        CvmdScal(nPitches, -1.0, t2);
        CvmdMul(nPitches, t2, t4, Gamma + Mp);
        CvmdMul(nPitches, t1, t4, Gamma + Mp + nPitches);
    }
    else{
        /* lambda = a-sum(ROld.*phi,1); */
        CvmdInit(nPitches*(iOrder-2), 0.0, t4);
        for(k = 0 ; k < iOrder-2 ; ++k){
            CvmdAdd(nPitches, tf[iOrder-2-k], hf[k+iOrder-2], t1);
            CvmdMul(nPitches, phi + k*nPitchesOld, t1, t2);
            CvmdAddInplace(nPitches, t4, t2);
        }
        CvmdAdd(nPitches, tf[iOrder-1], hf[iOrder-3], t2); //a is now in t2
        CvmdSub(nPitches, t2, t4, lambda);
        
        /* mu = -sum(ROld.*psi,1); */
        CvmdInit(nPitches*(iOrder-2), 0.0, mu);
        for(k = 0 ; k < iOrder-2 ; ++k){
            CvmdAdd(nPitches, tf[iOrder-2-k], hf[k+iOrder-2], t1);
            CvmdMul(nPitches, psi + k*nPitchesOld, t1, t2);
            CvmdSubInplace(nPitches, mu, t2);
        }
        
        
        /* phi = [phi;zeros(1,nPitches)]+(ones(iOrder-1,1)*lambda).*gammaNew; */
        CvmdCopy((iOrder-2)*nPitchesOld, phi, t1); //could copy less
        for(k = 0 ; k < iOrder-2 ; ++k){
            CvmdMul(nPitches, lambda, Gamma + Mp*(iOrder-2) + nPitchesOld*k, t2);
            CvmdAdd(nPitches, t1 + nPitchesOld*k, t2, phi + nPitches*k);
        }
        CvmdMul(nPitches, lambda,
                Gamma + Mp*(iOrder-2) + nPitchesOld*(iOrder-2),
                phi + nPitches*(iOrder-2));
        
        /* psi = [psi;zeros(1,nPitches)]+(ones(iOrder-1,1)*mu).*gammaNew; */
        CvmdCopy((iOrder-2)*nPitchesOld, psi, t1); //could copy less
        for(k = 0 ; k < iOrder-2 ; ++k){
            CvmdMul(nPitches, mu, Gamma + Mp*(iOrder-2) + nPitchesOld*k, t2);
            CvmdAdd(nPitches, t1 + nPitchesOld*k, t2, psi + nPitches*k);
        }
        CvmdMul(nPitches, mu,
                Gamma + Mp*(iOrder-2) + nPitchesOld*(iOrder-2),
                psi + nPitches*(iOrder-2));
        
        
        /* alphaNew = sum(RNew(1:end-1,:).*gammaNew,1); */
        CvmdInit(nPitches, 0.0, t1);
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdAdd(nPitches, tf[iOrder-1-k], hf[k + iOrder-1], t2);
            CvmdMul(nPitches, t2, Gamma + Mp*(iOrder-2) + nPitchesOld*k, t3);
            CvmdAddInplace(nPitches, t1, t3);
        }
        
        /*     b = (ones(iOrder-1,1)*(alphaOld-alphaNew)).*gammaNew+...
         [zeros(1,nPitches);gammaNew(1:iOrder-2,:)]+...
         [gammaNew(2:end,:);zeros(1,nPitches)]-...
         [gammaOld(1:iOrder-2,:);zeros(1,nPitches)]+...
         (ones(iOrder-1,1)*psi(end,:)).*phi-...
         (ones(iOrder-1,1)*phi(end,:)).*psi;
         */
        CvmdSub(nPitches, alpha, t1, t2);
        //make alpha1 -> alphaOld in the next iteration
        CvmdCopy(nPitches, t1, alpha);
        
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdMul(nPitches, t2, Gamma + Mp*(iOrder-2) + nPitchesOld*k,
                    t3 + k*nPitches);
            if(k != 0){
                CvmdAddInplace(nPitches, t3 + k*nPitches,
                               Gamma + Mp*(iOrder-2) + nPitchesOld*(k-1));
            }
            if(k != iOrder - 2){
                CvmdAddInplace(nPitches, t3 + k*nPitches,
                               Gamma + Mp*(iOrder-2) + nPitchesOld*(k+1));
                
                CvmdSubInplace(nPitches, t3 + k*nPitches,
                               Gamma + Mp*(iOrder-3) + nPitchesOldOld*k);
            }
            CvmdMul(nPitches, psi + nPitches*(iOrder-2), phi + nPitches*k, t1);
            CvmdAddInplace(nPitches, t3 + k*nPitches, t1);
            CvmdMul(nPitches, phi + nPitches*(iOrder-2), psi + nPitches*k, t1);
            CvmdSubInplace(nPitches, t3 + k*nPitches, t1);
        }
        // b is now in t3
        
        /* nu = sum(RNew(1:end-1,:).*b)./gammaNew(end,:); */
        CvmdInit(nPitches, 0.0, t1);
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdAdd(nPitches, tf[iOrder-1-k], hf[k + iOrder-1], t2);
            CvmdMul(nPitches, t2, t3 + k*nPitches, t4);
            CvmdAddInplace(nPitches, t1, t4);
        }
        CvmdDiv(nPitches, t1, Gamma + Mp*(iOrder-2) + nPitchesOld*(iOrder-2), t2);
        //nu is now in t2
        
        /*      gammaOld = gammaNew;
         gammaNew = nan(iOrder,nPitches);
         gammaNew(iOrder,:) = 1./(nu+RNew(iOrder,:));
         gammaNew(1:iOrder-1,:) = (ones(iOrder-1,1)*...
         (gammaNew(iOrder,:)./gammaOld(end,:))).*b;
         */
        CvmdAdd(nPitches, tf[0], hf[2*(iOrder-1)], t1);
        CvmdAddInplace(nPitches, t2, t1);
        CvmdInverse(nPitches, t2, Gamma + Mp*(iOrder-1) + nPitches*(iOrder-1));
        CvmdDiv(nPitches, Gamma + Mp*(iOrder-1) + nPitches*(iOrder-1),
                Gamma + Mp*(iOrder-2) + nPitchesOld*(iOrder-2), t4);
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdMul(nPitches, t3 + nPitches*k, t4,
                    Gamma + Mp*(iOrder-1) + nPitches*k);
        }
    }
}

/*
 Implements a special version for the system
 
 R = T - H
 
 where
 
 t_2 = h_0, t_3 = h_1, ....
 
 In this case phi, psi, mu, lambda and a can be left out.
 */

void single_pitch::update_gamma_p(int iOrder, int nPitches, int nPitchesOld, 
                                  int nPitchesOldOld,  FTYPE * alpha,
                                  FTYPE * Gamma){
    
    if( iOrder == 1 ){
        CvmdSub(nPitches, tf[0], hf[0], t2);
        CvmdInverse(nPitches, t2, Gamma);
    }
    else if( iOrder == 2 ){
        CvmdSub(nPitches, tf[0], hf[0], t1);     // R00 = t[0] - h[0];
        CvmdSub(nPitches, tf[1], hf[1], t2); // R01 = t[1] - h[1];
        CvmdSub(nPitches, tf[0], hf[2], t3); // R11 = t[0] - h[2];
        CvmdMul(nPitches, t1, t3, t3);
        CvmdMul(nPitches, t2, t2, t4);
        CvmdSubInplace(nPitches, t3, t4);
        CvmdInverse(nPitches, t3, t4);
        
        CvmdMul(nPitches, t2, Gamma, alpha);
        
        CvmdScal(nPitches, -1.0, t2);
        CvmdMul(nPitches, t2, t4, Gamma + Mp);
        CvmdMul(nPitches, t1, t4, Gamma + Mp + nPitches);
    }
    else{
        
        /* alphaNew = sum(RNew(1:end-1,:).*gammaNew,1); */
        CvmdInit(nPitches, 0.0, t1);
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdSub(nPitches, tf[iOrder-1-k], hf[k + iOrder-1], t2);
            CvmdMul(nPitches, t2, Gamma + Mp*(iOrder-2) + nPitchesOld*k, t3);
            CvmdAddInplace(nPitches, t1, t3);
        }
        
        /*     b = (ones(iOrder-1,1)*(alphaOld-alphaNew)).*gammaNew+...
         [zeros(1,nPitches);gammaNew(1:iOrder-2,:)]+...
         [gammaNew(2:end,:);zeros(1,nPitches)]-...
         [gammaOld(1:iOrder-2,:);zeros(1,nPitches)]+...
         */
        CvmdSub(nPitches, alpha, t1, t2);
        //make alpha1 -> alphaOld in the next iteration
        CvmdCopy(nPitches, t1, alpha);
        
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdMul(nPitches, t2, Gamma + Mp*(iOrder-2) + nPitchesOld*k,
                    t3 + k*nPitches);
            if(k != 0){
                CvmdAddInplace(nPitches, t3 + k*nPitches,
                               Gamma + Mp*(iOrder-2) + nPitchesOld*(k-1));
            }
            if(k != iOrder - 2){
                CvmdAddInplace(nPitches, t3 + k*nPitches,
                               Gamma + Mp*(iOrder-2) + nPitchesOld*(k+1));
                
                CvmdSubInplace(nPitches, t3 + k*nPitches,
                               Gamma + Mp*(iOrder-3) + nPitchesOldOld*k);
            }
        }
        // b is now in t3
        
        /* nu = sum(RNew(1:end-1,:).*b)./gammaNew(end,:); */
        CvmdInit(nPitches, 0.0, t1);
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdSub(nPitches, tf[iOrder-1-k], hf[k + iOrder-1], t2);
            CvmdMul(nPitches, t2, t3 + k*nPitches, t4);
            CvmdAddInplace(nPitches, t1, t4);
        }
        CvmdDiv(nPitches, t1, Gamma + Mp*(iOrder-2) + nPitchesOld*(iOrder-2), t2);
        //nu is now in t2
        
        /*      gammaOld = gammaNew;
         gammaNew = nan(iOrder,nPitches);
         gammaNew(iOrder,:) = 1./(nu+RNew(iOrder,:));
         gammaNew(1:iOrder-1,:) = (ones(iOrder-1,1)*...
         (gammaNew(iOrder,:)./gammaOld(end,:))).*b;
         */
        CvmdSub(nPitches, tf[0], hf[2*(iOrder-1)], t1);
        CvmdAddInplace(nPitches, t2, t1);
        CvmdInverse(nPitches, t2, Gamma + Mp*(iOrder-1) + nPitches*(iOrder-1));
        CvmdDiv(nPitches, Gamma + Mp*(iOrder-1) + nPitches*(iOrder-1),
                Gamma + Mp*(iOrder-2) + nPitchesOld*(iOrder-2), t4);
        for(k = 0 ; k < iOrder-1 ; ++k){
            CvmdMul(nPitches, t3 + nPitches*k, t4,
                    Gamma + Mp*(iOrder-1) + nPitches*k);
        }
    }
}


FTYPE single_pitch::compute_obj(FTYPE omega, FTYPE * x, 
                                int iOrder, FTYPE * ac, FTYPE * as){
    
    FTYPE * ccVectors = vector(2*(iOrder+1));
    
    /*
     crossCorrelationVectors = 0.5*[nData;...
     sin((1:2*iOrder+1)'*0.5*omega*nData)./...
     (sin((1:2*iOrder+1)'*0.5*omega))];
     */
    ccVectors[0] = nData*(1+T0_REG);
    
#ifdef MKL
    CvmdSinRange(2*iOrder+1, range+1, 0.5*omega*nData, t4, t1);
    CvmdSinRange(2*iOrder+1, range+1, 0.5*omega, t4, t2);
#else
    CvmdSinRange(1, 2*iOrder+1, 0.5*omega*nData, t1);
    CvmdSinRange(1, 2*iOrder+1, 0.5*omega, t2);
#endif
    CvmdDiv(2*iOrder+1, t1, t2, ccVectors + 1);
    CvmdScal(2*(iOrder+1), 0.5, ccVectors);
    
    /* n = -(nData-1)/2:-(nData-1)/2+nData-1; */
    FTYPE * n = vector(nData);
    CvmdAddConstant(nData, -(nData-1)*0.5, range, n);
    
    /* C = cos(omega*n'*(1:iOrder));
     S = sin(omega*n'*(1:iOrder));
     
     bc = C'*x;
     bs = S'*x;
     */
    FTYPE * bc = vector(iOrder);
    FTYPE * bs = vector(iOrder);
    
#ifdef CHEBYSHEV
    /* Implementation using Chebyshev recursive method */
    FTYPE * cn = vector(nData);
    FTYPE * sn = vector(nData);
    FTYPE * cnm1 = vector(nData);
    FTYPE * snm1 = vector(nData);
    FTYPE * cnm2 = vector(nData);
    FTYPE * snm2 = vector(nData);
    FTYPE * tN = vector(nData);
    FTYPE * tw;
    
    CvmdCosSinRange(nData, n, omega, tN, cn, sn);
    bc[0] = CvmdDot(nData, cn, x);
    bs[0] = CvmdDot(nData, sn, x);
    
    if( iOrder > 1 ){
        CvmdInit(nData, 1.0, cnm2);
        CvmdMul(nData, cn, cn, tN);
        CvmdAxpy(nData, -2, tN, cnm2);
        CvmdScal(nData, -1.0, cnm2);
        bc[1] = CvmdDot(nData, cnm2, x);
        
        CvmdMul(nData, cn, sn, snm2);
        CvmdScal(nData, 2, snm2);
        bs[1] = CvmdDot(nData, snm2, x);
        
        if( iOrder > 2 ){
            CvmdCopy(nData, cnm2, cnm1);
            CvmdCopy(nData, snm2, snm1);
            CvmdCopy(nData, cn, cnm2);
            CvmdCopy(nData, sn, snm2);
            for( int l = 3 ; l <= iOrder ; ++l ){
                CvmdMul(nData, cn, cnm1, tN);
                CvmdAxpy(nData, -2, tN, cnm2);
                CvmdScal(nData, -1, cnm2);
                bc[l-1] = CvmdDot(nData, cnm2, x);
                
                CvmdMul(nData, cn, snm1, tN);
                CvmdAxpy(nData, -2, tN, snm2);
                CvmdScal(nData, -1, snm2);
                bs[l-1] = CvmdDot(nData, snm2, x);
                
                /* update using swap of pointers */
                tw = cnm1;
                cnm1 = cnm2;
                cnm2 = tw;
                
                tw = snm1;
                snm1 = snm2;
                snm2 = tw;
            }
        }
    }
    del_vector(cn);
    del_vector(sn);
    del_vector(cnm1);
    del_vector(snm1);
    del_vector(cnm2);
    del_vector(snm2);
    del_vector(tN);
#else
    // standard impl
    for( int l = 1 ; l <= iOrder ; ++l ){
        CvmdCosSinRange(nData, n, omega*l, t1, t2, t3);
        bc[l-1] = CvmdDot(nData, t2, x);
        bs[l-1] = CvmdDot(nData, t3, x);
    }
#endif
    
    /*
     t = crossCorrelationVectors(1:iOrder);
     h = crossCorrelationVectors(3:2*iOrder+1);
     */
    
    FTYPE * t = ccVectors;
    FTYPE * h = vector(2*(iOrder-1)+1);
    CvmdCopy(2*(iOrder-1)+1, ccVectors + 2, h);
    
    /*
     T = toeplitz(t);
     H = hankel(h(1:iOrder), h(iOrder:end));
     
     ac = (T + H)\bc;
     as = (T - H)\bs;
     */
    
    /* Solve the system using the TH algorithm */
    int K = ((iOrder+1)*(iOrder))/2;
    FTYPE * gamma = vector(K);
    
    th(iOrder-1, t, h, gamma);
    solve(iOrder-1, t, h, bc, gamma, ac);
    
    for(int k = 0; k < 2*(iOrder-1)+1 ; ++k)
        h[k] = -h[k];
    
    thp(iOrder-1, t, h, gamma);
    solve(iOrder-1, t, h, bs, gamma, as);
    
    /* obj = bc'*ac + bs'*as */
    FTYPE obj = CvmdDot(iOrder, bc, ac) + CvmdDot(iOrder, bs, as);
    
    del_vector(ccVectors);
    del_vector(n);
    del_vector(bc);
    del_vector(bs);
    del_vector(h);
    del_vector(gamma);
    
    return obj;
}


/*
 
 Compute the objective at omega, for data x, with length nData
 and order iOrder
 */
FTYPE single_pitch::compute_obj(FTYPE omega, FTYPE * x, int iOrder){
    
    FTYPE * ac = new FTYPE [iOrder];
    FTYPE * as = new FTYPE [iOrder];
    
    FTYPE obj = compute_obj(omega, x, iOrder, ac, as);
    
    delete [] ac;
    delete [] as;
    
    return obj;
}


/* returns the maximum objective at order iOrder for the 
 last computed objective
 */
FTYPE single_pitch::max_obj(int iOrder){
    
    FTYPE max =  -std::numeric_limits<float>::max();
    FTYPE * obj = costFunctionMatrix[iOrder-1];
    
    for( int k = 0 ; k < nPitchesAll[iOrder-1] ; ++k){
        max = MAX(obj[k], max);
    }
    
    return max;
}


/* returns the omega that maximes the objective at order iOrder
 for the lates computed objetive
 */
FTYPE single_pitch::argmax_obj(int iOrder){
    //      modified by Liming Shi
    FTYPE max_value = -std::numeric_limits<float>::max();
    int arg;
    FTYPE * obj = costFunctionMatrix[iOrder-1];
    
    for( int k = 0 ; k < nPitchesAll[iOrder-1] ; ++k){
        if( obj[k] > max_value){
            max_value = obj[k];
            arg =  k;
        }
    }
    return (2*M_PI*(arg+minFftIndex))/nFftGrid;
}


/* refine of the grid solution that was latest computed 
 returns a pointer to an internal length maxModelOrder array
 
 eps is the size of of the last interval measured in radian per samples
 */
FTYPE * single_pitch::refine(FTYPE *x, FTYPE eps){
    FTYPE res = 2*M_PI/nFftGrid;
    FTYPE omega, omega_l, omega_u;
    /* eps accuracy could be modified to follow CRLB or % of a note */
    
    for(int iOrder = 1 ; iOrder <= maxModelOrder ; ++iOrder){
        omega = argmax_obj(iOrder);
        omega_l = omega - res;
        omega_u = omega + res;
        omega_0h[iOrder-1] = golden(x, iOrder, omega_l, omega_u, eps);
    }
    return omega_0h;
}

void single_pitch::compute_max_on_grid(void){
    
    for(int iOrder = 1 ; iOrder <= maxModelOrder ; ++iOrder){
        omega_0h[iOrder-1] = argmax_obj(iOrder);
    }
}

/* refine of the grid solution that was latest computed 
 returns an estimate in radians per sample
 
 eps is the size of of the last interval measured in radian per samples
 */
FTYPE single_pitch::refine_single(FTYPE *x, int iOrder, FTYPE eps){
    FTYPE res = 2*M_PI/nFftGrid;
    FTYPE omega, omega_l, omega_u;
    /* eps accuracy could be modified to follow CRLB or % of a note */
    
    if (iOrder >= 1){
        omega = argmax_obj(iOrder);
        omega_l = omega - res;
        omega_u = omega + res;
        return golden(x, iOrder, omega_l, omega_u, eps);
    }
    return -1;
}

/* refine of the grid solution that was latest computed 
 returns a pointer to an internal length maxModelOrder array
 
 default eps = 1e-4
 */
FTYPE * single_pitch::refine(FTYPE *x){
    return refine(x, 1e-4);
}

/* 
 Applies the golden section search to the objective defined
 via the function compute_obj()
 
 The folden section search follows the description in
 
 A. Antoniou and W.-S. Lu
 Practical Optimization
 Algorithms and Engineering Applications
 2007
 Springer
 */
FTYPE single_pitch::golden(FTYPE * x, int iOrder, 
                           FTYPE omega_L, FTYPE omega_U, FTYPE eps){
    FTYPE K = 1.618033988749895;
    
    FTYPE omega_l = omega_L;
    FTYPE omega_u = omega_U;
    
    FTYPE Ik = (omega_U - omega_L)/K;
    
    FTYPE omega_a = omega_u - Ik;
    FTYPE omega_b = omega_l + Ik;
    
    FTYPE fa = -compute_obj(omega_a, x, iOrder);
    FTYPE fb = -compute_obj(omega_b, x, iOrder);
    
    while( Ik > eps or omega_u < omega_l ){
        
        Ik = Ik/K;
        
        if( fa >= fb ){
            omega_l = omega_a;
            omega_a = omega_b;
            omega_b = omega_l + Ik;
            fa = fb;
            fb = -compute_obj(omega_b, x, iOrder);
        }
        else{
            omega_u = omega_b;
            omega_b = omega_a;
            omega_a = omega_u - Ik;
            fb = fa;
            fa = -compute_obj(omega_a, x, iOrder);
        }
    }
    
    if( fa > fb )
        return 0.5*(omega_b + omega_u);
    else if( fa < fb )
        return 0.5*(omega_l + omega_a);
    else
        return 0.5*(omega_a + omega_b);
}

/* perform model order selection on the data x at 
 the frequencies currently located in omega_h
 as obtained using the refine procedure
 
 This model order selection method is outlined in
 
 J.K. Nielsen, M.G. Christensen, A.T. Cemgil and S.H. Jensen
 Bayesian Model Comparision with the g-Prior
 IEEE Transaction on Signal Processing
 Vol 62, pp. 225-238, 2014
 
 */
// int single_pitch::model_order_selection(FTYPE * x){
//   return single_pitch::model_order_selection(x, 0);
// }



int single_pitch::model_order_selection(){
    
    int order;
    int maxFftIndex = (int)floor(nFftGrid*pitchBounds[1]);
    int minPitchGrid = minFftIndex;
    int maxPitchGrid = maxFftIndex;
    int nPitches = maxFftIndex - minFftIndex + 1;
    //     double temp;
    
    
    double max_temp =  -std::numeric_limits<float>::max();
    float deltaFreq= ( (float) maxFftIndex-minFftIndex)/(nPitches-1)/nFftGrid;
    FTYPE * obj;
    FTYPE * logMarginalLikelihood= vector(maxModelOrder);
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        logMarginalLikelihood[ii-1]=0.0;
        obj = costFunctionMatrix[ii-1];
        for( int k = 0 ; k < M ; ++k){
            if( obj[k] > max_temp){
                max_temp = obj[k];
            }
        }
        for( int jj = 1 ; jj <= M ; ++jj){
            //             temp=;
            logMarginalLikelihood[ii-1]= logMarginalLikelihood[ii-1]+exp(costFunctionMatrix[ii-1][jj-1]-max_temp);
        }
        logMarginalLikelihood[ii-1]=log(logMarginalLikelihood[ii-1])+max_temp+log(deltaFreq);
        max_temp =  -std::numeric_limits<float>::max();
    }
    FTYPE * joint_logMarginalLikelihood= vector(maxModelOrder+1);
    joint_logMarginalLikelihood[0]=0.0;
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        joint_logMarginalLikelihood[ii]=logMarginalLikelihood[ii-1];
    }
    del_vector(logMarginalLikelihood);
    //     compute the posterior PMF;
    for( int ii = 1 ; ii <= maxModelOrder+1 ; ++ii){
        joint_logMarginalLikelihood[ii-1]=joint_logMarginalLikelihood[ii-1]+logModelPrior[ii-1];
    }
    
    int arg;
    max_temp= -std::numeric_limits<float>::max();
    obj = joint_logMarginalLikelihood;
    for( int k = 0 ; k < maxModelOrder+1; ++k){
        if( obj[k] > max_temp){
            max_temp = obj[k];
            arg=k;
        }
    }
    order=arg;
    del_vector(joint_logMarginalLikelihood);
    return order;
}
/*
 The single step function that performs the three steps
 
 1. Calculate the objective function for all candidate model
 order and on the Fourier grid
 2. Refine the best estimates for each model order
 3. Perform model order selection
 
 and return the estimated frequency in radians per sample
 */
FTYPE single_pitch::est(FTYPE * x){
    
    compute(x);
    refine(x, 1e-6);
    estOrder = model_order_selection();
    
    if( estOrder >= 0 )
        return omega_0h[estOrder-1];
    else
        return 0.0;
}

FTYPE single_pitch::est(FTYPE * x, FTYPE lnBFZeroOrder, FTYPE eps){
    
    compute(x);
    refine(x, eps);
    estOrder = model_order_selection();
    
    if( estOrder >= 0 )
        return omega_0h[estOrder-1];
    else
        return 0.0;
}

/*
 The single step function that performs the three steps
 
 1. Calculate the objective function for all candidate model
 order and on the Fourier grid
 2. Perform model order selection
 3. Refine the for the selected model order
 
 and return the estimated frequency in radians per sample
 */
FTYPE * single_pitch::est_fast(FTYPE * x){
    
    
    compute(x);
    FTYPE x_square=0.0;
    for (int i = 1; i <= nData ; i++){
        x_square=x_square+pow(x[i-1],2.0);
    }
    for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
        for( int jj = 1 ; jj <= M ; ++jj){
            costFunctionMatrix[ii-1][jj-1]=costFunctionMatrix[ii-1][jj-1]*(1.0/(x_square+5e-3));
        }
    }
    //
    FTYPE v=1.0;
    double log_scale;
//    double unvoicing_scaled_alpha;
    FTYPE delta=3;
    FTYPE null_modellike=1.0;
    FTYPE w;
    FTYPE unvoicing_bar_alpha;
    FTYPE u=FTYPE (nData)/2.0;
    FTYPE * a = vector(M);
    FTYPE * b = vector(M);
    FTYPE * gHat = vector(M);
    FTYPE * tauVar = vector(M);
    FTYPE temp_sum=0.0;
    for( int iOrder = 1 ; iOrder <= maxModelOrder ; ++iOrder){
        // compute the Laplace Parameters
        w=(nData-2*iOrder-delta)/2.0;
        for( int jj = 1 ; jj <= M ; ++jj){
            a[jj-1]=(1-costFunctionMatrix[iOrder-1][jj-1])*(v+w-u);
            b[jj-1]=(u-v)*costFunctionMatrix[iOrder-1][jj-1]+2*v+w-u;
            gHat[jj-1]=(b[jj-1]+sqrt(b[jj-1]*b[jj-1]-4*a[jj-1]*v))/(-2*a[jj-1]);
            tauVar[jj-1] = 1.0/(gHat[jj-1]*(1-costFunctionMatrix[iOrder-1][jj-1])*u/pow(1.0+gHat[jj-1]*(1-costFunctionMatrix[iOrder-1][jj-1]),2.0)-gHat[jj-1]*w/pow(1+gHat[jj-1],2.0));
            //         compute the pitch loglikelihood
            costFunctionMatrix[iOrder-1][jj-1]=log(gHat[jj-1]*(delta-2)/2.0)+(nData-2*iOrder-delta)/2.0*log(1+gHat[jj-1])-nData/2.0*log(1+gHat[jj-1]*(1-costFunctionMatrix[iOrder-1][jj-1]))+1.0/2.0*log(2.0*M_PI*tauVar[jj-1]);
        }
    }
    del_vector(a);
    del_vector(b);
    del_vector(gHat);
    del_vector(tauVar);
    
    //     obtain the unnormalized posterior
    if (isnan(scaled_alpha_buffer[0][0])){
        for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
            for( int jj = 1 ; jj <= M ; ++jj){
                bar_alpha[ii-1][jj-1]=logpi[jj-1]+(costFunctionMatrix[ii-1][jj-1]);
            }
        }
        unvoicing_bar_alpha=log(logModelPrior[1]*null_modellike);
    }else{
        

        //        obtain matrix multiplication obj.B'*obj.scaled_alpha_buffer
        for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
            for( int jj = 1 ; jj <= M ; ++jj){
                temp_sum=0.0;
                for ( int kk = 1 ; kk <= maxModelOrder; ++kk){
                    temp_sum=temp_sum+B[kk-1][ii-1]*(scaled_alpha_buffer[kk-1][jj-1]);
                }
                bar_alpha[ii-1][jj-1]=temp_sum;
            }
        }

        
        
        
        
        //       obj.C(0,0)*obj.B'*obj.scaled_alpha_buffer*obj.A
//        using scaled_alpha_buffer for saving the data, since it will not be used
        for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
            for( int jj = 1 ; jj <= M ; ++jj){
                temp_sum=0.0;
                for ( int kk = 1 ; kk <= M; ++kk){
                    temp_sum=temp_sum+bar_alpha[ii-1][kk-1]*A[kk-1][jj-1];
                }
                scaled_alpha_buffer[ii-1][jj-1]=C[0][0]*temp_sum;
            }
        }
        
        for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
            for( int jj = 1 ; jj <= M ; ++jj){
                bar_alpha[ii-1][jj-1]= scaled_alpha_buffer[ii-1][jj-1];
            }
        }
        
        
        for( int ii = 1 ; ii <= maxModelOrder ; ++ii){
            for( int jj = 1 ; jj <= M ; ++jj){
                bar_alpha[ii-1][jj-1]=bar_alpha[ii-1][jj-1]+scaled_alpha_buffer2[ii-1][jj-1]*C[1][0]*unvoicing_scaled_alpha_buffer;
                bar_alpha[ii-1][jj-1]=log(bar_alpha[ii-1][jj-1])+(costFunctionMatrix[ii-1][jj-1]);
            }
        }

        unvoicing_bar_alpha=log(C[0][1]*(1-unvoicing_scaled_alpha_buffer)+unvoicing_scaled_alpha_buffer*C[1][1]);

    }
    
    
    
//    normalize the posterior
    log_scale= log_sumsum_exp( bar_alpha,maxModelOrder,M,unvoicing_bar_alpha);
    unvoicing_scaled_alpha=norm_prob( bar_alpha,maxModelOrder,M,unvoicing_bar_alpha,log_scale);
    

    int estOrder;
    int estimatedPitch_inx;
    FTYPE estimatedPitch;
    FTYPE max_value = -std::numeric_limits<float>::max();
    
    
    for (int ii = 0 ; ii < maxModelOrder ; ii++){
        for (int jj = 0 ; jj < M ; jj++){
            if (bar_alpha[ii][jj] > max_value){
                max_value = bar_alpha[ii][jj];
                estOrder= ii+1;
                estimatedPitch_inx=jj;
                
            }
            scaled_alpha_buffer[ii][jj]= exp(bar_alpha[ii][jj]);
        }
    }
    if (unvoicing_scaled_alpha < log(0.5))
    {
        log_scale= log(1-exp(unvoicing_scaled_alpha));
        norm_prob2( bar_alpha,maxModelOrder,M,log_scale,scaled_alpha_buffer2);
    }
    estimatedPitch=(2*M_PI*(estimatedPitch_inx+minFftIndex))/nFftGrid;
    unvoicing_scaled_alpha_buffer=exp(unvoicing_scaled_alpha);
    
    returned_value[0]=estimatedPitch;
    returned_value[1]=estOrder;
    returned_value[2]=1-unvoicing_scaled_alpha_buffer;
    return returned_value;

    
    
}


//  added by Liming Shi
inline double single_pitch::squared(double x)
{
    return x*x;
}

double single_pitch::normpdf(double x, double u, double s) { 
    return (1.0/(s*sqrt(2.0*M_PI)))*exp(-0.5*squared(x-u)/squared(s));
}

double single_pitch::log_sum_exp(double arr[], int count) 
{
    double maxVal = arr[0];
    double sum = 0;
    
    for (int i = 1 ; i < count ; i++){
        if (arr[i] > maxVal){
            maxVal = arr[i];
        }
    }
    
    for (int i = 0; i < count ; i++){
        sum += exp(arr[i] - maxVal);
    }
    return log(sum) + maxVal;
    
}

double single_pitch::norm_prob(FTYPE **arr, int row_num, int col_num, double A_num,double scale_par)
{
    for (int ii = 0 ; ii < row_num ; ii++){
        for (int jj = 0 ; jj < col_num ; jj++){
//            if (arr[ii][jj] > maxVal){
                arr[ii][jj]=arr[ii][jj]-scale_par;
//            }
        }
    }
//    if (A_num > maxVal){
        A_num = A_num-scale_par;
//    }
    

    return A_num;
    
}

int single_pitch::norm_prob2(FTYPE **arr, int row_num, int col_num, double scale_par,FTYPE **arrtarget)
{
    for (int ii = 0 ; ii < row_num ; ii++){
        for (int jj = 0 ; jj < col_num ; jj++){
            //            if (arr[ii][jj] > maxVal){
            arrtarget[ii][jj]=exp(arr[ii][jj]-scale_par);
            //            }
        }
    }
    //    if (A_num > maxVal){
//    A_num = A_num-scale_par;
    //    }
    
    
    return 0;
    
}


double single_pitch::log_sumsum_exp(FTYPE **arr, int row_num, int col_num, double A_num)
{
    double maxVal = arr[0][0];
    double sum = 0;
    for (int ii = 0 ; ii < row_num ; ii++){
        for (int jj = 0 ; jj < col_num ; jj++){
            if (arr[ii][jj] > maxVal){
                maxVal = arr[ii][jj];
            }
        }
    }
    if (A_num > maxVal){
        maxVal = A_num;
    }
    
    for (int ii = 0 ; ii < row_num ; ii++){
        for (int jj = 0 ; jj < col_num ; jj++){
            sum += exp(arr[ii][jj] - maxVal);
        }
    }
    sum += exp(A_num - maxVal);
    return log(sum) + maxVal;
    
}



double single_pitch::log_sum(double arr[], int count) 
{
    double sum = 0;
    
    for (int i = 0; i < count ; i++){
        sum += arr[i];
    }
    return log(sum);
    
}
double single_pitch::sum_vec(double arr[], int count) 
{
    double sum = 0;
    
    for (int i = 0; i < count ; i++){
        sum += arr[i];
    }
    return sum;
}

//double ramp(double x)
//{
//    if (x > 0)
//        return x;
//    else
//        return 0;
//}

