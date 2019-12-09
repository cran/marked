#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_INTEGER(m);                                       // number of capture occasions
    DATA_MATRIX(ch);                                       // capture history matrix
    DATA_VECTOR(freq);                                     // capture history frequency; >1 for fixed effect only
    DATA_IVECTOR(frst);                                    // occasion first seen for each history
    DATA_IVECTOR(lst);                                     // occasion last seen for each history
    DATA_IVECTOR(loc);                                     // 0 or 1, 1 if lost on capture at last event
    DATA_MATRIX(tint);                                     // time interval between occasions 

    DATA_MATRIX(phi_fixedDM);                              // Phi DM
    DATA_INTEGER(phi_nre);                                 // number of random effects for phi
    DATA_INTEGER(phi_krand);                               // number of columns in phi random effect DM
    DATA_MATRIX(phi_randDM);                               // phi random effect DM
    DATA_IMATRIX(phi_randIndex);                           // phi random effect indices for DM
    DATA_IVECTOR(phi_counts);                              // count of phi random effect indices by id
    DATA_IMATRIX(phi_idIndex);                             // phi random effect indices by id

    DATA_MATRIX(p_fixedDM);                               // p DM
    DATA_INTEGER(p_nre);                                  // number of random effects for p
    DATA_INTEGER(p_krand);                                // number of columns in p random effect DM
    DATA_MATRIX(p_randDM);                                // p random effect DM
    DATA_IMATRIX(p_randIndex);                            // p random effect indices for DM; index into phi_u
    DATA_IVECTOR(p_counts);                               // count of p random effect indices by id
    DATA_IMATRIX(p_idIndex);                              // p random effect indices by id; index into u_phi to construct phi_u

    DATA_INTEGER(getreals);                               // =1 if real phi's and p's are to be computed
    
    DATA_INTEGER(prior);                                  // either 0 or 1 (if specified prior means/sigma for betas and sigmas for random effects)
    DATA_VECTOR(mu_phi_prior);
    DATA_VECTOR(sigma_phi_prior);
    DATA_VECTOR(mu_p_prior);
    DATA_VECTOR(sigma_p_prior);
    DATA_VECTOR(random_mu_phi_prior);
    DATA_VECTOR(random_sigma_phi_prior);
    DATA_VECTOR(random_mu_p_prior);
    DATA_VECTOR(random_sigma_p_prior);
       
    PARAMETER_VECTOR(phi_beta);
    PARAMETER_VECTOR(p_beta);
    PARAMETER_VECTOR(log_sigma_phi);
    PARAMETER_VECTOR(log_sigma_p);
    PARAMETER_VECTOR(u_phi);
    PARAMETER_VECTOR(u_p);
    
    Type f=0;
    int n=ch.rows();                  // number of capture histories
    int kphi=phi_fixedDM.cols()-1;    // number of Phi DM columns
    int kp=p_fixedDM.cols()-1;        // number of p DM columns

    if(phi_krand>0)	                                        // likelihood contribution for n(0,1) re for 
       for (int i=0;i<=phi_nre-1;i++)	   	     
           f-= dnorm(u_phi(i),Type(0),Type(1),true);
           
    if(p_krand>0)	                                        // likelihood contribution for n(0,1) re for p
        for (int i=0;i<=p_nre-1;i++)
	      f-= dnorm(u_p(i),Type(0),Type(1),true);
	      
    if(prior>0) 
    {  
    	for(int i=0;i<=kphi-1;i++) 
       		f-= dnorm(phi_beta(i), mu_phi_prior(i), sigma_phi_prior(i), true);
    	for(int i=0;i<=kp-1;i++) 
       		f-= dnorm(p_beta(i), mu_p_prior(i), sigma_p_prior(i), true);
       	if(phi_krand>0)	
	    	for(int i=0;i<=phi_krand-1;i++) 
    	   		f-= dnorm(exp(log_sigma_phi(i)), random_mu_phi_prior(i), random_sigma_phi_prior(i), true) +log_sigma_phi(i);
       	if(p_krand>0)	
    		for(int i=0;i<=p_krand-1;i++) 
       			f-= dnorm(exp(log_sigma_p(i)), random_mu_p_prior(i), random_sigma_p_prior(i), true)+ log_sigma_p(i);
    }
    int nphicounts=n;                                       // number of counts for phi random effects by id
    if(phi_nre==0)nphicounts=0;
    int npcounts=n;                                         // number of counts for p random effects by id
    if(p_nre==0)npcounts=0;
  
    for(int i=1;i<=n;i++)                // loop over capture histories - one per animal
    {
       vector<Type> phi(m-1);            // temp vector for Phis for each occasion for a single history
       vector<Type> p(m-1);              // temp vector for Phis for each occasion for a single history
       vector<Type> phicumprod(m);       // cummulative survival probability across occasions
       vector<Type> cump(m);             // cummulative probability of being seen across occasions
       Type pch;                         // probability of capture history
       int i1,i2,j,L;	                 // miscellaneous ints
       Type mu;                          // link function value
	   
       vector<Type> p_u(p_idIndex.cols());        // define random effects vector for p and Phi used                  
       vector<Type> phi_u(phi_idIndex.cols());    // just for this capture history copied from full vectors
       p_u.setZero();	   
       phi_u.setZero();	   
	     phi.setZero();
	     p.setZero();
 	     if(nphicounts >0)                          // if any random effects for phi, copy values from u_phi to phi_u
	     {
   		    if(phi_counts(i-1)==0)
		         phi_u(0)=0;
		      else
		         for(j=0;j<=phi_counts(i-1)-1;j++)
			          phi_u(j)=u_phi(phi_idIndex(i-1,j)-1);
	     } 

       if(npcounts >0)                           // if any random effects for p, copy values from u_phi to phi_u
	     {
   		    if(p_counts(i-1)==0)
		         p_u(0)=0;
		      else
		         for(j=0;j<=p_counts(i-1)-1;j++)
		            p_u(j)=u_p(p_idIndex(i-1,j)-1);
	      } 
		
   	    for(j=0;j<=m-1;j++)              // loop over all occasions              
        {   
 	          phicumprod(j)=1.0;             // set cummulative survival to 1
            cump(j)=1.0;                   // set cummulative p to 1
        }	
       
        i1=(m-1)*(i-1);                   // compute beginning index in design matrix
       
	      for(j=frst(i-1)+1;j<=m;j++)       // loop over occasions from frst to m
	      {
	        i2=i1+j-1;                      // increment index in design matrix

         /////// phi computation //////////
	        if(phi_fixedDM(i2-1,kphi)== -1)
	        {
	           mu=0;
	           for(L=1;L<=kphi;L++)
                 mu+=phi_fixedDM(i2-1,L-1)*phi_beta(L-1);    // fixed portion of mean 
			       if(nphicounts>0)
			          if(phi_counts(i-1) > 0)	                        // random portion of mean if any
	              {
	                 for(L=1;L<=phi_krand;L++)
		                  if(phi_randIndex(i2-1,L-1)>0)
	                       mu+=phi_randDM(i2-1,L-1)*phi_u(phi_randIndex(i2-1,L-1)-1)*exp(log_sigma_phi(L-1));
                 }	
              phi(j-2)=1/(1+exp(-mu));                        // compute phi for the interval
              if(tint(i-1,j-2)!=1)
                phi(j-2)=pow(phi(j-2),tint(i-1,j-2));       // adjust phi for the time interval length
	          } else
	              phi(j-2)=phi_fixedDM(i2-1,kphi);                // real fixed value for this phi
	        
 		  /////// p computation /////////
	      if(p_fixedDM(i2-1,kp)== -1)
	      {
   	        mu=0;
	          for(L=1;L<=kp;L++)
                mu+=p_fixedDM(i2-1,L-1)*p_beta(L-1);         // fixed portion of mean 
			      if(npcounts>0)
			         if(p_counts(i-1) > 0)	                         // random portion of mean if any
	             {
	                for(L=1;L<=p_krand;L++)
                     if(p_randIndex(i2-1,L-1)>0)
	                       mu+=p_randDM(i2-1,L-1)*p_u(p_randIndex(i2-1,L-1)-1)*exp(log_sigma_p(L-1));
                }	
                p(j-2)=1/(1+exp(-mu));                                     
          } else
	           p(j-2)=p_fixedDM(i2-1,kp);                      // real fixed value for this p
			 
	       phicumprod(j-1)=phicumprod(j-2)*phi(j-2);                             // compute cummulative survival
	       cump(j-1)=cump(j-2)*((1-p(j-2))*(1-ch(i-1,j-1))+p(j-2)*ch(i-1,j-1));  // compute cummulative capture probability
       } 
       pch=0.0;                                                      // initialize prob of the capture history to 0
       for(j=lst(i-1);j<=((1-loc(i-1))*m+loc(i-1)*lst(i-1));j++)     // loop over last occasion to m unless loss on capture
       {                                                             // to compute prob of the capture history 
          if((loc(i-1)==1) | (j==m))
             pch=pch+cump(j-1)*phicumprod(j-1);                      // probability of history given possible last occasion alive
          else       
             pch=pch+cump(j-1)*phicumprod(j-1)*(1-phi(j-1));         // probability of history given possible last occasion alive
       }   
       f-= log(pch+1E-24)*freq(i-1);                                 // sum log-likelihood log(pr(ch))
       
       if(getreals==1)                                               // if requested report phi and p real values
	     {
	        ADREPORT(phi);
	        ADREPORT(p);
	     }    
	}
    return f;
}
