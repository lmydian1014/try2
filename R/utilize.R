### given the data, compute all the intermediate value


get_nii_grids<- function(img){
  return(list(X = 1:img@dim_[2]*img@srow_x[1]+img@srow_x[4],
              Y = 1:img@dim_[3]*img@srow_y[2]+img@srow_y[4],
              Z = 1:img@dim_[4]*img@srow_z[3]+img@srow_z[4]))
}

get_nii_coords<- function(img){
  grids = get_nii_grids(img)
  coords = expand.grid(grids)
  return(coords)
}


classify.accuracy = function(est_type,true_type,levels=c(1,0)){
  if(is.null(levels))
    levels = unique(true_type)
  tab = table(factor(est_type,levels=levels), factor(true_type,levels=levels))
  TP = tab[1,1]
  FP = tab[1,2]
  TN = tab[2,2]
  FN = tab[2,1]
  Sensitivity = TP/(TP+FN)
  Specificity = TN/(TN+FP)
  PPV = TP/(TP+FP)
  NPV = TN/(TN+FN)
  Accuracy = (TP+TN)/sum(tab)
  FDR = 1 - PPV
  return(c(TP=TP,FP=FP,TN=TN,FN=FN,
           Sensitivity=Sensitivity,
           Specificity=Specificity,
           PPV=PPV,NPV=NPV,FDR=FDR,
           Accuracy=Accuracy))
}

reg_single_v = function(v, modality, Z, Y_1, Y_2, X){
  Yv = c()
  Xv = c()
  if(modality == 1){
    for(i in c(1:n)){
      Yv = c(Yv, Y_1[v,i])
      Xv = c(Xv, X[v,i])
    }
  }
  else{
    for(i in c(1:n)){
      Yv = c(Yv, Y_2[v,i])
      Xv = c(Xv, X[v,i])
    }
  }
  M = cbind(Xv, Z)
  fit = lm(Yv ~ M)
  return(c(fit$residual))
}



G_thres = function(x, thres){
	return(ifelse(x > thres, x, 0))
}

################################ loglikelihood
loglikelihood = function(n, e_pos, e_neg, c_star, tau_1_sq, tau_2_sq, Xmat, Y_pos, Y_neg, thres){
	Sc = Xmat %*% c_star
	Se_pos = Xmat %*% e_pos
	Se_neg = Xmat %*% e_neg
	mu_pos = c(G_thres(Sc, thres)) * Se_pos
	mu_neg = c(G_thres(-Sc, thres)) * Se_neg
	r = (tau_1_sq - tau_2_sq)/(tau_1_sq + tau_2_sq) 	## vector
	s_sq = (tau_1_sq + tau_2_sq)/4	
	H1 = rowSums((Y_pos - mu_pos)^2 + (Y_neg - mu_neg)^2 - 2 * c(r) * (Y_pos - mu_pos) * (Y_neg - mu_neg))
	H2 = 2 * (1 - r^2) * s_sq
	H3 = 2 * pi * s_sq * sqrt(1 - r^2)
	logL = -n * log(H3) - H1/H2
	return(sum(logL))
}

loglike = function(n, e_pos, e_neg, c_star, tau_1_sq, tau_2_sq, Xmat, Y_1, Y_2, thres, rho){
  tau_1 = sqrt(tau_1_sq)
  tau_2 = sqrt(tau_2_sq)
  Sc = Xmat %*% c_star
  Se_pos = Xmat %*% e_pos
  Se_neg = Xmat %*% e_neg
  eta_pos = c(G_thres(Sc, thres)) * Se_pos
  eta_neg = c(G_thres(-Sc, thres)) * Se_neg
  mu1 = eta_pos + eta_neg; mu2 = eta_pos - eta_neg
  H1 = 2 * pi * tau_1 * tau_2 * sqrt(1-rho^2)
  H2 = rowSums((Y_1 - mu1)/tau_1) - 2 * rowSums(rho * ((Y_1 - mu1)/tau_1) * ((Y_2 - mu2)/tau_2)) + rowSums((Y_2 - mu2)/tau_2)
  H3 = -n * log(H1) - H2/(2*(1-rho^2))
  return(sum(H3))
}

LLH_design = function(dat){
  tau_1 = c(sqrt(dat$tau_1_sq)); tau_2 = c(sqrt(dat$tau_2_sq)); 
  eta_pos = dat$eta_pos; eta_neg = dat$eta_neg;
  Y1 = dat$Y_1;  Y2 = dat$Y_2;  
  rho = c(dat$rho); 
  mu1 = eta_pos + eta_neg; mu2 = eta_pos - eta_neg
  H1 = 2 * pi * tau_1 * tau_2 * sqrt(1-rho^2)
  H2 = rowSums((Y1 - mu1)/tau_1) - 2 * rowSums(rho * ((Y1 - mu1)/tau_1) * ((Y2 - mu2)/tau_2)) + rowSums((Y2 - mu2)/tau_2)
  H3 = -n * log(H1) - H2/(2*(1-rho^2))
  return(sum(H3))
}

spectral_kernel = function(d, a1 = 1/2, a2 = 1/4, a3 = 1/4, l1 = 0.1, l2 = 0.1, l3 = 0.2, mu1 = 10, mu2 = 10, mu3 = 10){
    a1 * exp(- d^2 / l1^2) * cos(mu1*d) + a2 * exp(- d^2 / l2^2) * cos(mu2*d) + a3 * exp(- d^2 / l3^2) * cos(mu3*d)
}



approx_eigen = function(n, grids, num, l = 0.05){

  D = sqrt(distance(grids))
  K = (1 + sqrt(3) * D/l) * exp(- sqrt(3) * D/l)
 #K = spectral_kernel(D, a1 = 1/8, a2 = 3/8, a3 = 4/8, l1 = 0.1, l2 = 0.1, l3 = 0.2, mu1 = 50, mu2 = 10, mu3 = 10)
  #K = spectral_kernel(D, a1 = 1/8, a2 = 3/8, a3 = 4/8, l1 = 0.1, l2 =0.1, l3 = 0.2, mu1 = 0.5, mu2 = 0.10, mu3 = 0.1)
  #K = exp(- D / l)
  V = nrow(K)
  eig = eigen(K)
  #eig = eigs_sym(K, num)
  val_K = eig$values
  vec_K = eig$vectors
  lambda = val_K/V
  print(dim(K))
  print(dim(vec_K))
  Xmat = t(t(K %*% vec_K) * sqrt(n)/val_K)
  var_ratio = sum(lambda[c(1:num)])/sum(lambda)
  return(list('var_ratio' = var_ratio, 'lambda' = lambda[c(1:num)], 'Xmat' = Xmat[, c(1:num)]))
}

#### figure plot
onefig.levelplot = function(z,x,y,titles="fig1",cols=blue2red_cols){
  pdat = list(z = z,group=rep(titles,each=length(z)),x=rep(x,times=1),y=rep(y,times=1))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1)
  
}



onefig.levelplot.points = function(z,x,y,pts,titles="fig1",cols=blue2red_cols,pch=19,
                                   pcol="black",pcex=1.5){
  pdat = list(z = z,group=rep(titles,each=length(z)),x=rep(x,times=1),y=rep(y,times=1))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1,
            panel=function(...){
              panel.levelplot(...)
              lpoints(pts[,1], pts[,2], 
                      pch=pch,col=pcol,cex=pcex)
            })
  
}



twofigs.levelplot = function(z1,z2,x,y,titles=c("fig1","fig2"),cols=blue2red_cols){
  pdat = list(z = c(z1,z2),group=rep(titles,each=length(z1)),x=rep(x,times=2),y=rep(y,times=2))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1)
}


twofigs.levelplot.lines = function(z1,z2,x,y,lx,ly,lwd=2,lcol="black",titles=c("fig1","fig2"),cols=blue2red_cols){
  pdat = list(z = c(z1,z2),group=rep(titles,each=length(z1)),x=rep(x,times=2),y=rep(y,times=2))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1,
            panel=function(...){
              panel.levelplot(...)
              llines(lx,ly, lwd=lwd,col=lcol)
            })
}


# xgrid = expand.grid(seq(0,1,length=50),seq(0,1,length=100))
# z1 = exp(-10*(xgrid[,1]-xgrid[,2])^2)-0.8
# z2 = 1.5*exp(-100*(xgrid[,1]-xgrid[,2])^2)-0.8
# lx = seq(0.5,1,length=100)
# ly = lx^2
#twofigs.levelplot.lines(z1,z2,x=xgrid[,1],y=xgrid[,2],lx=lx,ly=ly)


threefigs.levelplot = function(z1,z2,z3,x,y,titles=c("fig1","fig2","fig3"),cols=blue2red_cols){
  pdat = list(z = c(z1,z2,z3),group=rep(titles,each=length(z1)),x=rep(x,times=3),y=rep(y,times=3))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1)
}


fourfigs.levelplot = function(z1,z2,z3,z4,x,y,titles=c("fig1","fig2","fig3","fig4"),cols=blue2red_cols,
                              layout = NULL,xlab=NULL,ylab=NULL,panel=panel.levelplot){
  pdat = list(z = c(z1,z2,z3,z4),group=factor(rep(titles,each=length(z1)),levels=titles),x=rep(x,times=4),y=rep(y,times=4))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1,
            layout=layout,xlab=xlab,ylab=ylab,panel=panel)
}


multi.figs.levelplot = function(zmat,x,y,titles = NULL,cols=blue2red_cols, layout = NULL,xlab=NULL,ylab=NULL,panel=panel.levelplot,
                                colorkey=TRUE,alpha.regions=1,at=NULL){
  if(is.null(titles)){
    titles = paste("fig",1:ncol(zmat))
  }
  pdat = list(z = c(zmat),group = factor(rep(titles,each=nrow(zmat)),levels=titles),x=rep(x,times=ncol(zmat)),y=rep(y,times=ncol(zmat)))
  if(is.null(at)){
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1,
            layout=layout,xlab=xlab,
            ylab=ylab,panel=panel,alpha.regions=alpha.regions,colorkey=colorkey)
  }
  else{
    levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1,
              layout=layout,xlab=xlab,
              ylab=ylab,panel=panel,alpha.regions=alpha.regions,colorkey=colorkey,at=at)
  }
}


blue2red_cols = c("#000080", "#000083", "#000087", "#00008B", "#00008F", "#000093", "#000097", "#00009B",
            "#00009F", "#0000A3", "#0000A7", "#0000AB", "#0000AF", "#0000B3", "#0000B7", "#0000BB",
            "#0000BF", "#0000C3", "#0000C7", "#0000CB", "#0000CF", "#0000D3", "#0000D7", "#0000DB",
            "#0000DF", "#0000E3", "#0000E7", "#0000EB", "#0000EF", "#0000F3", "#0000F7", "#0000FB",
            "#0004FF", "#0008FF", "#000CFF", "#0010FF", "#0014FF", "#0018FF", "#001CFF", "#0020FF",
            "#0024FF", "#0028FF", "#002CFF", "#0030FF", "#0034FF", "#0038FF", "#003CFF", "#0040FF",
            "#0044FF", "#0048FF", "#004CFF", "#0050FF", "#0054FF", "#0058FF", "#005CFF", "#0060FF",
            "#0064FF", "#0068FF", "#006CFF", "#0070FF", "#0074FF", "#0078FF", "#007CFF", "#0080FF",
            "#0083FF", "#0087FF", "#008BFF", "#008FFF", "#0093FF", "#0097FF", "#009BFF", "#009FFF",
            "#00A3FF", "#00A7FF", "#00ABFF", "#00AFFF", "#00B3FF", "#00B7FF", "#00BBFF", "#00BFFF",
            "#00C3FF", "#00C7FF", "#00CBFF", "#00CFFF", "#00D3FF", "#00D7FF", "#00DBFF", "#00DFFF",
            "#00E3FF", "#00E7FF", "#00EBFF", "#00EFFF", "#00F3FF", "#00F7FF", "#00FBFF", "#00FFFF",
            "#04FFFB", "#08FFF7", "#0CFFF3", "#10FFEF", "#14FFEB", "#18FFE7", "#1CFFE3", "#20FFDF",
            "#24FFDB", "#28FFD7", "#2CFFD3", "#30FFCF", "#34FFCB", "#38FFC7", "#3CFFC3", "#40FFBF",
            "#44FFBB", "#48FFB7", "#4CFFB3", "#50FFAF", "#54FFAB", "#58FFA7", "#5CFFA3", "#60FF9F",
            "#64FF9B", "#68FF97", "#6CFF93", "#70FF8F", "#74FF8B", "#78FF87", "#7CFF83", "#80FF80",
            "#83FF7C", "#87FF78", "#8BFF74", "#8FFF70", "#93FF6C", "#97FF68", "#9BFF64", "#9FFF60",
            "#A3FF5C", "#A7FF58", "#ABFF54", "#AFFF50", "#B3FF4C", "#B7FF48", "#BBFF44", "#BFFF40",
            "#C3FF3C", "#C7FF38", "#CBFF34", "#CFFF30", "#D3FF2C", "#D7FF28", "#DBFF24", "#DFFF20",
            "#E3FF1C", "#E7FF18", "#EBFF14", "#EFFF10", "#F3FF0C", "#F7FF08", "#FBFF04", "#FFFF00",
            "#FFFB00", "#FFF700", "#FFF300", "#FFEF00", "#FFEB00", "#FFE700", "#FFE300", "#FFDF00",
            "#FFDB00", "#FFD700", "#FFD300", "#FFCF00", "#FFCB00", "#FFC700", "#FFC300", "#FFBF00",
            "#FFBB00", "#FFB700", "#FFB300", "#FFAF00", "#FFAB00", "#FFA700", "#FFA300", "#FF9F00",
            "#FF9B00", "#FF9700", "#FF9300", "#FF8F00", "#FF8B00", "#FF8700", "#FF8300", "#FF8000",
            "#FF7C00", "#FF7800", "#FF7400", "#FF7000", "#FF6C00", "#FF6800", "#FF6400", "#FF6000",
            "#FF5C00", "#FF5800", "#FF5400", "#FF5000", "#FF4C00", "#FF4800", "#FF4400", "#FF4000",
            "#FF3C00", "#FF3800", "#FF3400", "#FF3000", "#FF2C00", "#FF2800", "#FF2400", "#FF2000",
            "#FF1C00", "#FF1800", "#FF1400", "#FF1000", "#FF0C00", "#FF0800", "#FF0400", "#FF0000",
            "#FB0000", "#F70000", "#F30000", "#EF0000", "#EB0000", "#E70000", "#E30000", "#DF0000",
            "#DB0000", "#D70000", "#D30000", "#CF0000", "#CB0000", "#C70000", "#C30000", "#BF0000",
            "#BB0000", "#B70000", "#B30000", "#AF0000", "#AB0000", "#A70000", "#A30000", "#9F0000",
            "#9B0000", "#970000", "#930000", "#8F0000", "#8B0000", "#870000", "#830000", "#800000")

