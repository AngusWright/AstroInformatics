library(shiny)
library(bslib)



vars <- setdiff(names(iris), "Species")

ui <- page_sidebar(
  sidebar = sidebar(open = "open", 
selectInput("model_type", "Model Type:",
            c("Linear" = "linear", 
              "Gaussian" = "gaussian",
              "Poissonian" = "poissonian",
              "Gamma" = "gamma"),
            selected='gaussian'),    
selectInput("noise_model", "Noise Model:",
            c("Constant" = "constant",
            "Linear" = "linear",
            "Log-Linear" = "log-linear"),
            selected='log-linear'),
selectInput("plot_type", "Figure Type:",
            c("Data & Model" = "data_model",
              "Model Bias" = "bias",
              "Bias-Width" = "correlation"),
            selected='data_model'),
checkboxInput("correct_htrsk", "Correct for Heteroskedasticity?", FALSE),
sliderInput("n.noise",
          "Number of noise realisations per model (logarithmic):",
          min = 0,
          max = 4,
          step=0.01,
          value = 2),
sliderInput("noise_sd",
            withMathJax("Noise range (logarithmic) :"),
            min = -3, 
            max = 1,
            step = 0.01, 
            value = c(-2,-1)),
sliderInput("model_sd",
            withMathJax("Range of model param (slope/stdev)"),
            min = 0.1, 
            max = 1.9,
            step = 0.01, 
            value = c(0.9,1.1)),
sliderInput("nmodels",
            expression("Number of models to simulate"),
            min = 2,
            max = 8,
            step = 1, 
            value = 3),
sliderInput("sample_range",
            expression("x-position sampling range"),
            min = -10,
            max = 10,
            step = 0.1, 
            value = c(-4,4)),
sliderInput("num_samples",
            expression("Number of x-position samples (logarithmic)"),
            min = 0,
            max = 4,
            step = 0.1, 
            value = 2.0),
sliderInput("annot_cex",
            expression("Annotation Scaling"),
            min = 0.5,
            max = 2.5,
            step = 0.1, 
            value = 1.5)
  ),
  plotOutput("plot", height="100%",width="100%", fill=TRUE)
)

server <- function(input, output, session) {
   selectedData <- reactive({
      iris[, c(input$xcol, input$ycol)]
    })
  
  clusters <- reactive({
    kmeans(selectedData(), input$clusters)
  })
  
  output$plot <- renderPlot({
    palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))
  
    model<-input$noise_model 
    min.noise.sd<-10^(input$noise_sd[1])
    max.noise.sd<-10^(input$noise_sd[2])
    #Define the grid of model widths
    sds<-seq(input$model_sd[1],input$model_sd[2],length=input$nmodels)
    if (input$model_type=='gaussian') { 
      #Define the truth model N(0,1)
      truth<-function(x) dnorm(x,mean=0,sd=1)
      #Define the fitted model dnorm N(mu,sig)
      fitted<-function(x,mean,width) dnorm(x,mean=mean,sd=width)
    } else if (input$model_type=='linear') {
      #Define the truth model y=x
      truth<-function(x) x
      #Define the fitted model y=mx+b
      fitted<-function(x,intercept,slope) slope*x+intercept
    } else { 
      stop(paste("model",input$model_type,"not yet implemented")) 
    }
    #Define the x sample points 
    x<-seq(input$sample_range[1],input$sample_range[2],length=10^input$num_samples)
    #Define the number of noise realisations
    n.noise<-10^input$n.noise

  amp<-fits<-matrix(NA,nrow=length(sds),ncol=n.noise)
  #Evaluate the truth at every x
  y<-truth(x)
  #Generate the noise 
  if (model=='constant') { 
    noise.model<-rep((min.noise.sd+max.noise.sd)/2,length(x))
  } else if (model=="linear") { 
    noise.model<-seq(min.noise.sd,max.noise.sd,len=length(x))
  } else if (model=="log-linear") { 
    noise.model<-10^seq(log10(min.noise.sd),log10(max.noise.sd),len=length(x))
  } else { 
    stop("unknown noise model") 
  }
  #Loop through model widths
  for (i in 1:length(sds)) {
  #Loop through noise realisations 
    for (j in 1:n.noise) {
  #Add noise to y data
      yobs<-y+rnorm(length(x),sd=noise.model)
      if (input$model_type=='gaussian') { 
         #Fit the model for amplitude and centroid
               tmp.fit<-optim(par=list(mu=-1,A=1), 
                   fn=function(par) return=sum((yobs-
                                                  par[2]*fitted(x,
                                                                mean=par[1],
                                                                width=sds[i])
                                                )^2/noise.model))
        #Save the best-fit parameters 
        fits[i,j]<-tmp.fit$par[1]
        amp[i,j]<-tmp.fit$par[2]
      } else if (input$model_type=='linear') { 
         #Fit the model for slope and intercept 
        tmp.fit<-optim(par=list(intercept=1), 
                   fn=function(par) return=sum((yobs-
                                                  fitted(x,
                                                         slope=sds[i],
                                                         intercept=par[1])
                                                )^2/noise.model))
        #Save the best-fit parameters 
        fits[i,j]<-tmp.fit$par[1]
      } 
      if (input$correct_htrsk) { 
        if (input$model_type=='gaussian') { 
          #Fit the model for amplitude and centroid with hetroskedasticity weights 
          resid<-abs(yobs-amp[i,j]*fitted(x,fits[i,j],sds[i]))
          resid_fit<-lm(resid~obs,
                        data=list(obs=fitted(x,mean=fits[i,j],width=sds[i]),
                                  resid=resid))
          resid_wgt<-1/predict(resid_fit)^2
          resid_wgt<-resid_wgt/sum(resid_wgt)
          tmp.fit<-optim(par=list(mu=-1,A=1), 
              fn=function(par) return=sum((yobs-par[2]*fitted(x,mean=par[1],width=sds[i]))^2/noise.model*resid_wgt))
          #Save the best-fit parameters 
          fits[i,j]<-tmp.fit$par[1]
          amp[i,j]<-tmp.fit$par[2]
        } else if (input$model_type=='linear') { 
          #Fit the model for amplitude and centroid with hetroskedasticity weights 
          resid<-abs(yobs-fitted(x,fits[i,j],sds[i]))
          resid_fit<-lm(resid~obs,
                        data=list(obs=fitted(x,intercept=fits[i,j],slope=sds[i]),
                                  resid=resid))
          resid_wgt<-1/predict(resid_fit)^2
          resid_wgt<-resid_wgt/sum(resid_wgt)
          tmp.fit<-optim(par=list(intercept=-1), 
              fn=function(par) return=sum((yobs-
                                             fitted(x,
                                                    intercept=par[1],
                                                    width=sds[i])
                                           )^2/noise.model*resid_wgt))
          #Save the best-fit parameters 
          fits[i,j]<-tmp.fit$par[1]
        }
      }
    }
  }
  #Initialise the plot
  par(mar = c(4.1, 4.1, 0, 1))
  if (input$plot_type=='data_model') { 
    if (input$model_type=='linear') { 
      magicaxis::magplot(x,yobs,pch=20,xlab='x',ylab='f(x)',ylim=range(x),
              cex.lab=input$annot_cex,cex.axis=input$annot_cex)
    } else { 
      magicaxis::magplot(x,yobs,pch=20,xlab='x',ylab='n(x)',ylim=c(-0.05,0.5),
              cex.lab=input$annot_cex,cex.axis=input$annot_cex)
    } 
    #Loop over model widths 
    for (i in 1:length(sds)) {
      #Plot 100 of the best-fit models 
      for (j in sample(x=1:n.noise,size=100)) { 
        if (input$model_type=='gaussian') { 
          lines(x,amp[i,j]*dnorm(x,mean=fits[i,j],sd=sds[i]),
                col=seqinr::col2alpha(RColorBrewer::brewer.pal(8,'Set2')[i],alpha=0.1),lwd=2)
        } else if (input$model_type=='linear') { 
          lines(x,sds[i]*x+fits[i,j],
                col=seqinr::col2alpha(RColorBrewer::brewer.pal(8,'Set2')[i],alpha=0.1),lwd=2)
        }
      }
    }
    #Add one noise realisation of the data
    points(x,yobs,pch=20)
    #Add the uncertanties
    magicaxis::magerr(x,yobs,ylo=noise.model)
    #Add the legend
    if (input$model_type=='linear') labbase<-'INT:'
    if (input$model_type=='gaussian') labbase<-'STD:'
    legend('topleft',lty=1,lwd=2,col=c('black',RColorBrewer::brewer.pal(8,'Set2')),
      legend=c(paste0(labbase,' 1.0'),paste(labbase,round(digits=2,sds[1:floor(length(sds)/2)]))),
      bty='n',pch=c(20,NA,NA),
      cex=input$annot_cex,inset=0.02)
    legend('topright',lty=1,lwd=2,
           col=c(RColorBrewer::brewer.pal(8,'Set2')[-(1:floor(length(sds)/2))]),
      legend=paste(labbase,round(digits=2,sds[-(1:floor(length(sds)/2))])),
                   bty='n',cex=input$annot_cex,inset=0.02)
    #Add the noise model
    #tcoord<-helpRfuncs::text.coord(loc='left',inset=0.15)
    tcoord<-c(quantile(x,prob=0.20),0.25)
    text(tcoord[1],tcoord[2],lab=paste0('Noise Model:\n',model),pos=3,cex=input$annot_cex) 
  } else if (input$plot_type=='bias') { 
    #Initialise the best-fit means plot
    magicaxis::magplot(density(fits[1,],bw=0.01/sqrt(12),kern='rect',
                    from=min(c(-0.2,min(fits))),
                    to=max(c(0.2,max(fits)))),
            xlab='Model Bias',ylab='PDF',type='n',
            cex.lab=input$annot_cex,cex.axis=input$annot_cex)
    #Loop over each of the model widths
    for (i in 1:length(sds)) {
      lines(density(fits[i,],bw=0.05/sqrt(12),kern='rect',
                    from=-2*max(abs(sds-1)),
                    to=2*max(abs(sds-1))),
            col=RColorBrewer::brewer.pal(8,'Set2')[i],lwd=2)
    }
  } else if (input$plot_type=='correlation') { 
    #Initialise the best-fit means plot
    if (input$model_type=='linear') labbase<-'Intercept'
    if (input$model_type=='gaussian') labbase<-'Stdev'
    magicaxis::magplot(jitter(rep(sds,n.noise)),fits,
            col=rep(RColorBrewer::brewer.pal(8,'Set2')[1:length(sds)],n.noise),
            ylab='Model Bias',xlab=paste('Model',labbase),type='p',pch=20,
            cex.lab=input$annot_cex,cex.axis=input$annot_cex)
    dat<-data.frame(sds=rep(sds,n.noise),fits=as.numeric(fits))
    new<-data.frame(sds=seq(min(sds)-1,max(sds)+1,len=1e3))
    lm.fit<-lm(fits~sds,dat)
    matplot(new$sds,
          cbind(predict(lm.fit,newdata=new,interval='confidence'),
                predict(lm.fit,newdata=new,interval='prediction')[,-1]),
          lwd=c(2,1,1,1,1),lty=c(1,2,2,3,3),col='black',type='l',add=TRUE)
  }

    #magicaxis::magplot(selectedData(),
    #     col = clusters()$cluster,
    #     pch = 20, cex = 3, side=1:4,labels=c(T,T,F,F), 
    #     xlab="xlabel", #gsub("\."," ",input$xcol),
    #     ylab="ylabel") #gsub("\."," ",input$ycol))
    #points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
  #} else { 
  #  par(mar = c(4.1, 4.1, 0, 1))
  #  magicaxis::magplot(selectedData(),
  #       col = clusters()$cluster,
  #       pch = 15, cex = 3, side=1:4,labels=c(T,T,F,F), 
  #       xlab="xlabel", #gsub("\."," ",input$xcol),
  #       ylab="ylabel") #gsub("\."," ",input$ycol))
  #  points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
  #}
  })
}


# Create Shiny app ----
shinyApp(ui = ui, server = server)
