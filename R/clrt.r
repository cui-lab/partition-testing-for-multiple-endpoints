#########################################################################################
###### This file provides functions implement the Consonant Likelihood Ratio Test  ######
###### for multiple endpoints problem.                                             ######
###### Copyright(c) 2011 Bushi Wang.                                               ######
###### Modified May 3, 2011, version 1.0.1                                         ######
###### This code is distributed WITHOUT ANY WARRANTY.                              ######
#########################################################################################
###### Modified May 24, 2011, version 1.0.2                                        ######
###### In this version, .CLRT.PS4.test has been modified to return p-value and cv. ######
#########################################################################################
###### Modified September 19, 2011, version 1.1.0                                  ######
###### Build R package clrt.                                                       ######
#########################################################################################


# The program is only for the case with one primary and one secondary endpoint problem.
# Comparing k-1 doses of study drug with a control.

# If m is the number of endpoints, we restrict m=2.

# The natural partition divides the parameter space into (m+1)^(k-1)=26 pieces for k=4. 
#   Thus there are 26 null hypotheses to test. The null distribution of the parameter
#   of interest in all 26 hypotheses are identical, which assumes there is no dose effect
#   for both primary and secondary endpoints.

# We first generate simu=10^5 datasets, each contains k=4 groups of 2-by-1 vectors,
#   each group represent the patients receiving a paticular dose. Thus the number
#   of vectors in each group is the number of patients in this group.

# Next we calculate the LRT for each of these simu=10^5 datasets. Each LRT is a
#   26-by-1 vector indicating test statistic for each null hypothesis.

# In order to inplement the power of Consonant likelihood ratio test, we generate
#   powersimu=10^5 datasets following some senario with dose effect and known
#   correlation between endpoints. Then put each one of this powersimu=10^5
#   dataset into the original simulated dataset and calculate the p-value according
#   to the complete consonant adjustment. Thus we get powersimu=10^5 p-values each
#   one a 26-by-1 vector.

# For each of the powersimu=10^5 26-by-1 p-values, find the 2-by-k inference matrix.
#   Average among the powersimu=10^5 inference matrix we get the power.

.Constraint<-function(label,k=4){
	# label is a length k-1 vector indicating the label of an intersection hypothesis.
	# It follows the definition as in Table 1 and 2.
    hnullp<-hnulls<-matrix(0,nrow=k-1,ncol=k)
    for(i in 1:(k-1)){
        if(label[i]==1){
            hnullp[i,c(k,i)]<-c(1,-1)     # \Theta_i^{e_i}
        }else if(label[i]==2){
            hnulls[i,c(k,i)]<-c(1,-1)
        }
    }
    return(list('hnullp'=t(hnullp),'hnulls'=t(hnulls)))
}


.naturalpartition<-function(k=4,m=2,balanced=F){
    #k is the number of treatment arms, including control.
    nhypo<-(m+1)^(k-1)-1
    npmtrx<-matrix(NA,nrow=nhypo+1,ncol=k-1)
    for(a in 1:(k-1))npmtrx[,a]<-as.numeric(sapply(1:(m+1),function(x)rep(x,(m+1)^(k-1-a))))
    npmtrx<-npmtrx[1:nhypo,]
    # npmtrx gives all the possible labelings of intersection hypotheses
	
    if(balanced){
		# In the case of balanced treatment arms, the computation of critical constants can be greatly
		# 	simplified by sharing the same values among certain intersection hypohteses.
		# For example, hypothesis H_121 and H_112 will use the same critical value.
        labmtrx<-t(apply(npmtrx,1,sort))
        diff.label<-rep(1,nhypo)
        unique.label<-list()
        for(i in 1:nhypo){
            if(diff.label[i]==0){
                unique.label[[i]]<-0;
                next;
            }
            unique.label[[i]]<-which(sapply(1:nhypo,function(x)all(labmtrx[x,]==labmtrx[i,])))
            diff.label[ unique.label[[i]] ]<-0
        }
        balanced<-unique.label
    }
	
    # Some general definition, specified below.
	Comp.Amat<-list()
    component.array<-list()
    row.ind<-list()
    active.array<-list()
    if(k>3)col.ind<-lapply(1:(k-2),function(x)matrix(nrow=choose(k-1,x+1),ncol=choose(x+1,x)))
    
	for(cross in 1:(k-1)){
        component.array[[cross]]<-array(m+1,dim=c(m^cross,k-1,choose(k-1,cross)))
		# component.array also gives the labeling of all intersection hypotheses.
		# However, we may locate certain intersection hypothesis through component.array
        matrix.comb<-combn(k-1,cross)
        active.array[[cross]]<-matrix(NA,nrow=m^cross,ncol=cross)
        for(a in 1:cross)active.array[[cross]][,a]<-as.numeric(sapply(1:m,function(x)rep(x,m^(cross-a))))
        # active.array[[cross]] gives the possible combination of endpoints in a "cross"-way intersection.
		
		for(i in 1:choose(k-1,cross))component.array[[cross]][,matrix.comb[,i],i]<-active.array[[cross]]
        # component.array[[cross]] is an array, which consists of "choose(k-1,cross)" matrices. 
		#	The number "cross" indicates how many original hypotheses involved in the intersection.
		#	Each matrix, for example, component.array[[cross]][,,1], has each row a label for one paticular
		#	intersection. The number of rows in each matrix is m^cross.
		Comp.Amat[[cross]]<-list()
        # Comp.Amat are the constraint matrix for the component hypothesis.
		# It provide guidance on how to compute the restricted MLEs and the likelihood ratio statistics.
		# Comp.Amat[[1]] contain the constraint matrix for all the elementary hypotheses
		# Comp.Amat[[2]] contain the constraint matrix for all 2 way intersection hypotheses
		# 	ect... up to k-1.
		for(i in 1:choose(k-1,cross)){
            Comp.Amat[[cross]][[i]]<-list()
            for(j in 1:m^cross) Comp.Amat[[cross]][[i]][[j]]<-.Constraint(label=component.array[[cross]][j,,i],k=k)
        }
        
		# The list row.ind and col.ind tells how to perform the consonant adjustment.
		# Read the instruction in function .rNull.dist for definitions on ComponentLRT first.
		# Since each intersection hypothesis is located by one cell in the matrix of ComponentLRT[[1]][,,i],
		#	ComponentLRT[[2]][,,i], and ComponentLRT[[3]][,,i] for the ith replication, we give the index of those
		#	intersection hypotheses that needs to be tested in consonance adjustment.
		
		# For example, a three way intersection hypothesis H_121, will be located in ComponentLRT[[3]].
		#	According to active.array[[3]], it is on the 3rd row of ComponentLRT[[3]]. And ComponentLRT[[3]][,,i] only
		#	has one column. So the index for H_121 is ComponentLRT[[3]][3,1,i].
		# To reject H_121, we need to also reject at least two of H_123, H_131 and H_321. They are associated with three
		#	indices in ComponentLRT.
		# First, since they are all two-way intersections, they all located in ComponentLRT[[2]].
		#	row.ind[[2]][3,]=c(2,1,3) gives the three row index of these hypotheses,
		#	col.ind[[2]][1,]=c(1,2,3) gives the three col index of these hypotheses.
		#	They are ComponentLRT[[2]][2,1,i], ComponentLRT[[2]][1,2,i] and ComponentLRT[[2]][3,3,i].
		if(cross>1){
            row.ind[[cross-1]]<-matrix(NA,nrow=m^cross,ncol=cross)
            combn.mat<-combn(cross,cross-1)
            for(a in 1:m^cross)for(b in 1:cross)
                row.ind[[cross-1]][a,b]<-which(sapply(1:m^(cross-1),function(x)all(active.array[[cross-1]][x,1:(cross-1),drop=F]==active.array[[cross]][a,combn.mat[,b]])))
            combn.mat<-combn(k-1,cross-1)
            for(a in 1:nrow(col.ind[[cross-1]]))
                col.ind[[cross-1]][a,]<-which(sapply(1:choose(k-1,cross-1),function(x)all(combn.mat[,x]%in%combn(k-1,cross)[,a])))
        }
    }
    return(list('matrix'=npmtrx,'Comp.Amat'=Comp.Amat,'component.array'=component.array,'row.ind'=row.ind,'col.ind'=col.ind,'balanced'=balanced))
}


# LR compute the likelihood ratio statistics for each intersection hypothesis.
library(quadprog)
.LR<-function(Amat,k=4,m=2,Ybar,sigma,n=rep(50,4)){
    # Ybar is a m-by-k matrix, the last column in Ybar is the control.
	# sigma is a length m vector, includes stdev of each endpoint.
	# Amat is the constraint matrix for certain dose-endpoint combination from Comp.Amat
    sigma<-sum(n-1)*sigma^2
    ynullp<-solve.QP(Dmat=diag(n),dvec=diag(n)%*%Ybar[1,],Amat=Amat$hnullp)$solution
    ynulls<-solve.QP(Dmat=diag(n),dvec=diag(n)%*%Ybar[2,],Amat=Amat$hnulls)$solution
    LR<-( sigma[1]*sigma[2]
         /(sigma[1]+sum(n*(Ybar[1,]-ynullp)^2))/(sigma[2]+sum(n*(Ybar[2,]-ynulls)^2)) )^(-sum(n)/2)
    # To the power of -sum(n)/2 so that reject the null if LR > critical value.
    return(LR)
}


# .rNull.dist generate samples from the null distribution and compute their test statistics.
library(mvtnorm)
.rNull.dist<-function(n=rep(50,4),simu=10^5,rho=0){
    # n is the sample size for each dose. It is a vector of length 4, the last elements is the control.
    # The number of doses in this example, including control, is 4.
    # The number of endpoints is fixed at 2, for this function.
    k<-length(n)
    balanced<-(length(unique(n[1:(k-1)]))==1)
    NP<-.naturalpartition(k=k,m=2,balanced=balanced)
    # Generate null distribution
    sigma<-matrix(nrow=simu,ncol=2)
    Ybar<-array(dim=c(2,k,simu))
    ComponentLRT<-list()
    for(cross in 1:(k-1)) ComponentLRT[[cross]]<-array(dim=c(2^cross,choose(k-1,cross),simu))
    # ComponentLRT[[cross]] is a 2^cross * choose(k-1,cross) * simu array, where simu is the replication.
	# In the case of 4 treatment groups and 2 endpoints,
	#	ComponentLRT[[1]] is 2 by choose(3,1)=3 matrix, each cell is a one-way intersection hypothesis with
	#		index of the form x33, 3x3 or 33x. The x can be 1 or 2 for two endpoints, thus two rows.
	#	ComponentLRT[[2]] is 4 by choose(3,2)=3 matrix, each cell indicates a two-way intersection hypothesis
	#		with index of the form xx3, x3x or 3xx. The xx can be 11, 12, 21 and 22 from two endpoints, thus 
	#		there are 4 rows.
	#	ComponentLRT[[3]] is 8 by choose(3,3)=1 matrix, each cell indicates a three-way intersection hypothesis
	#		with index of the form xxx. The xxx comes for two endpoints combinations.
	
	# Each intersection hypothesis fits into one cell of the matrix in ComponentLRT[[1]][,,i], ComponentLRT[[2]][,,i],
	#	and ComponentLRT[[3]][,,i] for the ith replication. There are in total 26 of them.
	
	# The following piece of code comes from Normal random variable generating function.
	# To same computation time, part of the repeating computation is done here only once.
    retval<-list()
    VC<-matrix(rho,nrow=2,ncol=2)+(1-rho)*diag(2)
    for(i in 1:k){
        altersigma<-VC/n[i]
        ev<-eigen(altersigma, symmetric = TRUE)
        retval[[i]] <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    }

    for(i in 1:simu){
            sigma[i,]<-rchisq(2,df=sum(n-1))/sum(n-1)
            # sigma[i,] is the pooled estimate of the standard deviation, it is a vector
            #   of length 2, corresponding to each endpoint.

            Ybar[,,i]<-sapply(1:k,function(x)matrix(rnorm(2), nrow = 1) %*% retval[[x]])
            # 2 by 4 matrix, the last column is the control.

            for(cross in 1:(k-1))
                #ComponentLRT[[cross]][,,i]<-sapply(1:choose(k-1,cross),function(j)sapply(1:2^cross,function(x) .LR(Amat=NP$Comp.Amat[[cross]][[j]][[x]],k=k,m=2,Ybar=Ybar[,,i],sigma=sigma[i,],n=n)))
                for(j in 1:choose(k-1,cross))
                    ComponentLRT[[cross]][,j,i]<-sapply(1:2^cross,function(x) .LR(Amat=NP$Comp.Amat[[cross]][[j]][[x]],k=k,m=2,Ybar=Ybar[,,i],sigma=sigma[i,],n=n))
			# Compute all the likelihood ratio statistics. For some reason, it is faster to use a for loop.
    }
    return(list('ComponentLRT'=ComponentLRT,'NP'=NP,'n'=n,'simu'=simu))
}

#t<-proc.time()[1]
#set.seed(10)
#PS4.null<-.rNull.dist(n=rep(50,4),simu=10^3)
#proc.time()[1]-t

# .rAlternative.dist generate samples from alternative/true distribution and compute likelihood ratio statistics.
# It returns several statistics for different uses.
#	ComponentLRT is used by consonant likelihood ratio test
#	tstat is used by Dunnett gatekeeping procedure
#	DP.object is used by Decision path partition testing
#	NP, scenario, powersimu and delta is carried out to identify the problem.
.rAlternative.dist<-function(n=rep(50,4),powersimu=10^5,rho=0.3,delta=c(0,0),scenario=matrix(0,nrow=2,ncol=3)){
    # delta is the difference between treatment and control that we want to detect.
	k<-length(n)
    balanced<-(length(unique(n[1:(k-1)]))==1)
    NP<-.naturalpartition(k=k,m=2,balanced=balanced)

    sigma<-matrix(nrow=powersimu,ncol=2)
    Ybar<-array(dim=c(2,k,powersimu))
    ComponentLRT<-list()
	tstat<-array(dim=c(2,k-1,powersimu))
    # tstat is the t statistics used in Dunnett Gatekeeping procedure.
	for(cross in 1:(k-1)) ComponentLRT[[cross]]<-array(dim=c(2^cross,choose(k-1,cross),powersimu))
    if(length(delta)!=2)stop("delta should be of the same length as the number of endpoints!")
    trtmean<-cbind(scenario,delta)

    retval<-list()
    VC<-matrix(rho,nrow=2,ncol=2)+(1-rho)*diag(2)
    for(i in 1:k){
        altersigma<-VC/n[i]
        ev<-eigen(altersigma, symmetric = TRUE)
        retval[[i]] <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    }

    object<-list()
    for(i in 1:powersimu){
            sigma[i,]<-rchisq(2,df=sum(n-1))/sum(n-1)
            # sigma[i,] is the pooled estimate of the standard deviation, it is a vector
            #   of length 2, corresponding to each endpoint.

            Ybar[,,i]<-sapply(1:k,function(x) sweep(matrix(rnorm(2), nrow = 1) %*% retval[[x]], 2, trtmean[,x], "+") )
            # 2 by 4 matrix, the last column is the control

            for(cross in 1:(k-1))
                for(j in 1:choose(k-1,cross))
                    ComponentLRT[[cross]][,j,i]<-sapply(1:2^cross,function(x) .LR(Amat=NP$Comp.Amat[[cross]][[j]][[x]],k=k,m=2,Ybar=Ybar[,,i],sigma=sigma[i,],n=n))

    }
    return(list('ComponentLRT'=ComponentLRT,'NP'=NP,'scenario'=scenario,'delta'=delta,'powersimu'=powersimu))
}

#t<-proc.time()[1]
#PS4.alter<-.rAlternative.dist(n=rep(50,4),powersimu=10^4,rho=0.3,delta=c(0,0),scenario=matrix(c(0,0.65,0.65,0.65,0,0.65),nrow=2,ncol=3,byrow=T))
#proc.time()[1]-t


# In Windows System
#if(!is.loaded("consonant.adj")) dyn.load('./consonant.adj.dll')
# In Linux System
dyn.load("src/clrt.so")
.CLRT.PS4.cv<-function(alpha=0.025,null){
    # The number of endpoints comes from null, which is fixed at 2.
	n<-null$n
    simu<-null$simu
    k<-length(n)
    balanced<-(length(unique(n[1:(k-1)]))==1)
    NP<-null$NP
    ComponentLRT.cv<-list()
    for(cross in 1:(k-1)) ComponentLRT.cv[[cross]]<-matrix(nrow=2^cross,ncol=choose(k-1,cross))

    if(balanced){
    ComponentLRT.cv[[1]][1,]<-ComponentLRT.cv[[1]][2,]<-quantile(null$ComponentLRT[[1]],probs=1-alpha)
    for(cross in 2:(k-1)){
        n.adj<-dim(NP$col.ind[[cross-1]])[2]
        # List of variables needed by .C function
		# NP$col.ind[[cross-1]] is a choose(k-1,cross) by cross matrix
        # NP$row.ind[[cross-1]] is a 2^cross by cross matrix
		#	they share the same number of columns cross, which is the number of hypotheses to test in
		#	the consonance adjustment. Denoted by n.adj (equal to cross).
        # null$ComponentLRT[[cross-1]] is a 2^(cross-1) * choose(k-1,cross-1) * simu array
        # null$ComponentLRT[[cross]] is a 2^cross * choose(k-1,cross) * simu array
        for(a in 1:2^cross){
            cv<-sapply(1:choose(k-1,cross), function(j)
                            sapply(1:n.adj,function(x)ComponentLRT.cv[[cross-1]][NP$row.ind[[cross-1]][a,x],NP$col.ind[[cross-1]][j,x]]))
            # cv is a n.adj (cross) by choose(k-1,cross) matrix
			
			# The following .C code is the consonance adjustment, it will replace the likelihood ratio statistics of
			#	those intersection hypotheses do not satisfy consonance with 0.
            null$ComponentLRT[[cross]][a,,]<-.C("component_adj",
              as.double(null$ComponentLRT[[cross-1]]),
              as.integer(2^(cross-1)),
              as.integer(choose(k-1,cross-1)),
              ComponentLRT_post=as.double(null$ComponentLRT[[cross]][a,,]),
              as.integer(choose(k-1,cross)),
              as.integer(NP$row.ind[[cross-1]][a,]),
              as.integer(NP$col.ind[[cross-1]]),
              as.double(cv),
              as.integer(n.adj),
              as.integer(simu),
              PACKAGE="clrt" 
              )$ComponentLRT_post
        }
		# For balanced case, the critical constants are the same among intersection hypotheses with the same number
		#	of treatments involved regardless of endpoints involved.
        for(unq in k:1){
            L<-which(sapply(1:2^cross,function(x)length(unique(NP$component.array[[cross]][x,,1]))==unq))
            if(length(L)==0)next
            ComponentLRT.cv[[cross]][L,]<-quantile(null$ComponentLRT[[cross]][L,,],probs=1-alpha)
        }
    }
    }

    else{
    for(j in 1:(k-1))ComponentLRT.cv[[1]][,j]<-quantile(null$ComponentLRT[[1]][,j,],probs=1-alpha)
    for(cross in 2:(k-1)){
        n.adj<-dim(NP$col.ind[[cross-1]])[2]
        for(a in 1:2^cross){
            cv<-sapply(1:choose(k-1,cross), function(j)
                sapply(1:n.adj,function(x)ComponentLRT.cv[[cross-1]][NP$row.ind[[cross-1]][a,x],NP$col.ind[[cross-1]][j,x]]))
            null$ComponentLRT[[cross]][a,,]<-.C("component_adj",
                  as.double(null$ComponentLRT[[cross-1]]),
                  as.integer(2^(cross-1)),
                  as.integer(choose(k-1,cross-1)),
                  ComponentLRT_post=as.double(null$ComponentLRT[[cross]][a,,]),
                  as.integer(choose(k-1,cross)),
                  as.integer(NP$row.ind[[cross-1]][a,]),
                  as.integer(NP$col.ind[[cross-1]]),
                  as.double(cv),
                  as.integer(n.adj),
                  as.integer(simu),
                  PACKAGE="clrt" 
                  )$ComponentLRT_post
        }
        for(unq in k:1){
            L<-which(sapply(1:2^cross,function(x)length(unique(NP$component.array[[cross]][x,,1]))==unq))
            if(length(L)==0)next
            ComponentLRT.cv[[cross]][L,]<-sapply(1:choose(k-1,cross),function(j)quantile(null$ComponentLRT[[cross]][L,j,],probs=1-alpha))
        }
    }
    }

    return(ComponentLRT.cv)
}

#a<-proc.time()[1]
#cv<-.CLRT.PS4.cv(alpha=0.025,PS4.null)
#(proc.time()[1]-a)
#dyn.unload('./consonant.adj.dll')

.CLRT.PS4.test<-function(alpha=0.05,cv=NULL,n=rep(50,4),ybar,sigma,alter=NULL,powersimu=NULL,null=NULL,rho=0,delta=c(0,0),scenario=NULL){
    # The number of endpoints in this function is fixed at 2.
	# This function can take existing observations generated by .rAlternative.dist before, specify "alter".
	#	The delta and scenario will be obtained from the list "alter".
	# This function can also generate samples using .rAlternative.dist, give the desired sample size by "powersimu".
	#	And specify the correlation between endpoints with "rho" and the scenario want to test.
	# Otherwise, this function is used to test a specific observation, give their mean matrix and variance vector
	#	by "ybar" and "sigma". ybar is a m-by-k matrix treatment means, the last column is the control. sigma is
	#	a vector of length m, the stdev of each endpoint.
	k=length(n)
    m=2
    balanced<-(length(unique(n[1:(k-1)]))==1)
	
	if(is.null(cv))cv=.CLRT.PS4.cv(alpha=alpha,null=null)
	
    if(is.null(powersimu)){
        NP<-.naturalpartition(k=k,m=2,balanced=balanced)
        if(is.null(alter)){
            powersimu=1
            alter<-list('ComponentLRT'=list())
            for(cross in 1:(k-1)) alter$ComponentLRT[[cross]]<-array(dim=c(2^cross,choose(k-1,cross),1))
            ybar[,k]<-ybar[,k]+delta
            for(cross in 1:(k-1))
                for(j in 1:choose(k-1,cross))
                    alter$ComponentLRT[[cross]][,j,1]<-sapply(1:2^cross,function(x) .LR(Amat=NP$Comp.Amat[[cross]][[j]][[x]],k=k,m=2,Ybar=ybar,sigma=sigma,n=n))
        }
        else{
            powersimu<-alter$powersimu
            scenario<-alter$scenario
            delta<-alter$delta
        }
    }else{
        alter<-.rAlternative.dist(n=n,powersimu=powersimu,rho=rho,delta=delta,scenario=scenario)
        NP<-alter$NP
    }

    for(cross in 2:(k-1)){
        n.adj<-dim(NP$col.ind[[cross-1]])[2]
        for(a in 1:2^cross){
            cc<-sapply(1:choose(k-1,cross), function(j)
                sapply(1:n.adj,function(x)cv[[cross-1]][NP$row.ind[[cross-1]][a,x],NP$col.ind[[cross-1]][j,x]]))
            # The following .C code is the consonance adjustment, it will replace the likelihood ratio statistics of
			#	those intersection hypotheses do not satisfy consonance with 0.
            alter$ComponentLRT[[cross]][a,,]<-.C("component_adj",
                  as.double(alter$ComponentLRT[[cross-1]]),
                  as.integer(2^(cross-1)),
                  as.integer(choose(k-1,cross-1)),
                  ComponentLRT_post=as.double(alter$ComponentLRT[[cross]][a,,]),
                  as.integer(choose(k-1,cross)),
                  as.integer(NP$row.ind[[cross-1]][a,]),
                  as.integer(NP$col.ind[[cross-1]]),
                  as.double(cc),
                  as.integer(n.adj),
                  as.integer(powersimu),
                  PACKAGE="clrt" 
                  )$ComponentLRT_post
        }
    }
	# After the consonance adjustment, compare the likelihood ratio statistics with critical constants.
    if(is.null(null)){
		rejection<-matrix(nrow=(m+1)^(k-1)-1,ncol=powersimu)
		for(cross in 1:(k-1)){
			for(i in 1:choose(k-1,cross)){
				for(j in 1:m^cross){
					for(h in 1:((m+1)^(k-1)-1)){
						if(all(NP$matrix[h,]==NP$component.array[[cross]][j,,i])){
							rejection[h,]<-sapply(1:powersimu,function(x)alter$ComponentLRT[[cross]][j,i,x]>cv[[cross]][j,i])
							break
						}
					}
				}
			}
		}
	}
	# If p-value is desired, then compare likelihood ratio statistics with the null distribution and calculate p-value.
	else{
		pvalue<-rejection<-matrix(nrow=(m+1)^(k-1)-1,ncol=powersimu)
		for(cross in 1:(k-1)){
			for(i in 1:choose(k-1,cross)){
				for(j in 1:m^cross){
					for(h in 1:((m+1)^(k-1)-1)){
						if(all(NP$matrix[h,]==NP$component.array[[cross]][j,,i])){
							pvalue[h,]<-sapply(1:powersimu,function(x)mean(null$ComponentLRT[[cross]][j,i,]>alter$ComponentLRT[[cross]][j,i,x]))
							rejection[h,]<-pvalue[h,]<alpha
							break
						}
					}
				}
			}
		}
	}
	
    # The rejection rule follows Section 4.1
	Inference<-matrix(0,nrow=m,ncol=k-1)
    FWER<-numeric(powersimu)
    for(edp in 1:m){
        for(dose in 1:(k-1)){
            if(is.null(scenario)){
                Inference[edp,dose]<-all(rejection[which(NP$matrix[,dose]<=edp),1]==1)
                FWER<-NA
            }
            else if(scenario[edp,dose]==0) for(i in 1:powersimu){
                if(FWER[i]==1) Inference[edp,dose]<-Inference[edp,dose]+all(rejection[which(NP$matrix[,dose]<=edp),i]==1)
                else Inference[edp,dose]<-Inference[edp,dose]+(FWER[i]<-all(rejection[which(NP$matrix[,dose]<=edp),i]==1))
            }
            else for(i in 1:powersimu){
                Inference[edp,dose]<-Inference[edp,dose]+all(rejection[which(NP$matrix[,dose]<=edp),i]==1)
            }
        }
    }
    
	if(!is.null(null))
		return(list('alpha'=alpha,'Power'=Inference/powersimu,'FWER'=sum(FWER)/powersimu,'Partition'=apply(rejection,1,mean),'cv'=cv,'Pvalue'=pvalue,'Hypotheses'=NP$matrix[,1:3],'ComponentLRT'=alter$ComponentLRT ))
	else
		return(list('alpha'=alpha,'Power'=Inference/powersimu,'FWER'=sum(FWER)/powersimu,'Partition'=apply(rejection,1,mean),'cv'=cv,'Pvalue'=rejection,'Hypotheses'=NP$matrix[,1:3],'ComponentLRT'=alter$ComponentLRT ))
}

#a<-proc.time()[1]
#.CLRT.PS4.test(cv=cv,n=rep(50,4),powersimu=10^3,rho=0.3,delta=c(0,0),scenario=matrix(0.65,nrow=2,ncol=3))
#(proc.time()[1]-a)
#.CLRT.PS4.test(cv=cv,n=rep(50,4),ybar=matrix(8:1,nrow=2,ncol=4),sigma=c(2,4),delta=c(2,1))
#.CLRT.PS4.test(cv=cv,n=rep(50,4),alter=PS4.alter)

#source('./main.r') # Appended main.r source code above
if(!is.loaded("consonant.adj"))
# In Windows System
# dyn.load('./consonant.adj.dll')
# In Linux System
dyn.load("src/clrt.so")

clrt<-function(alpha=NULL,cv=NULL,n=NULL,grpmean=NULL,grpstdev=NULL,distr.null=NULL,distr.null.size=NULL,distr.alter=NULL,distr.alter.size=NULL,rho=0,delta=c(0,0),scenario=NULL,savenull=FALSE,savealter=FALSE,showpvalues=FALSE,showLRT=FALSE){
	if(is.null(alpha)){
		if(savenull){
			if(is.null(n))stop("Please specify sample size in each treatment groups and the control group.")
			if(is.null(distr.null.size))stop("Please specify how many samples are required for the null distribution.")
			distr.null=.rNull.dist(n=n,simu=distr.null.size,rho=0)
		}
		if(savealter){
			if(is.null(n))stop("Please specify sample size in each treatment groups and the control group.")
			if(is.null(distr.alter.size))stop("Please specify how many samples are required for the alternative distribution.")
			if(dim(scenario)[1]!=2)stop("This program can only handle cases with two endpoints.")
			if(dim(scenario)[2]!=length(n)-1)stop("Please specify the effect size for each treatment arm in each endpoint.")
			distr.alter=.rAlternative.dist(n=n,powersimu=distr.alter.size,rho=rho,delta=delta,scenario=scenario)
		}
		if(savenull&savealter) return(list('distr.null'=distr.null,'distr.alter'=distr.alter))
		else if(savenull) return(list('distr.null'=distr.null))
		else if(savealter) return(list('distr.alter'=distr.alter))
		else stop("Please specify significance level alpha if you want to perform test or calculate critical constants.")
	}
	else if(!is.null(grpmean) | !is.null(distr.alter) | !is.null(distr.alter.size)){
		if(is.null(cv)&is.null(distr.null)){
			if(is.null(n))stop("Please specify sample size in each treatment groups and the control group.")
			if(is.null(distr.null.size))stop("Please specify how many samples are required for the null distribution.")
			distr.null=.rNull.dist(n=n,simu=distr.null.size,rho=0)
		}
		if(!is.null(distr.null))showpvalue=TRUE
		else if(showpvalues){
			if(is.null(n))stop("Please specify sample size in each treatment groups and the control group.")
			if(is.null(distr.null.size))stop("Please specify how many samples are required for the null distribution.")
			distr.null=.rNull.dist(n=n,simu=distr.null.size,rho=0)
			showpvalue=TRUE
		}		
		test.result<-.CLRT.PS4.test(alpha=alpha,cv=cv,n=n,ybar=grpmean,sigma=grpstdev,alter=distr.alter,powersimu=distr.alter.size,null=distr.null,rho=rho,delta=delta,scenario=scenario)
		if(!savenull)distr.null=NULL
		if(!savealter)distr.alter=NULL
		if(!showLRT)test.result$ComponentLRT=NULL
		if(showpvalue) return(list('alpha'=alpha,'Power.Decision'=test.result$Power,'FWER'=test.result$FWER,'Partition.power'=test.result$Partition,'cv'=test.result$cv,'Partition.pvalue'=test.result$Pvalue,'Hypotheses'=test.result$Hypotheses,'ComponentLRT'=test.result$ComponentLRT,'distr.null'=distr.null,'distr.alter'=distr.alter))
		else return(list('alpha'=alpha,'Power.Decision'=test.result$Power,'FWER'=test.result$FWER,'Partition.power'=test.result$Partition,'cv'=test.result$cv,'Partition.rejection'=test.result$Pvalue,'Hypotheses'=test.result$Hypotheses,'ComponentLRT'=test.result$ComponentLRT,'distr.null'=distr.null,'distr.alter'=distr.alter))
	}
	else{
		if(is.null(distr.null)){
			if(is.null(n))stop("Please specify sample size in each treatment groups and the control group.")
			if(is.null(distr.null.size))stop("Please specify how many samples are required for the null distribution.")
			distr.null=.rNull.dist(n=n,simu=distr.null.size,rho=0)
			cv<-.CLRT.PS4.cv(alpha=alpha,null=distr.null)
		}
		else cv<-.CLRT.PS4.cv(alpha=alpha,null=distr.null)
		if(savenull)return(list('alpha'=alpha,'cv'=cv,'distr.null'=distr.null))
		else return(list('alpha'=alpha,'cv'=cv))
	}
}
