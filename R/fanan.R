#' This function to calculate the kruskal test(with neutrosophic data)
#'
#' @param dt ia a data frame
#'
#' @return kruskal test
#' @export
#' @importFrom stats kruskal.test
#' @importFrom utils capture.output
#' @examples
#'fac=c(rep("1",6),rep("2",6),rep("3",4))
#'t=c(0.4,0.42,0.04,0.46,0.08,0.33,0.13,0.003,0.0095,0.44,0.003,0.62,0.15,0.498,0.36,0.464)
#'i=c(0.06,0.071,0.5,0.14,0.03,0.30,0.45,0.074,0.17,0.28,0.48,0.072,0.62,0.148,0.831,0.761)
#'f=c(0.46,0.37,0.21,0.31,0.171,0.21,0.39,0.083,0.41,0.42,0.31,0.18,0.29,0.748,0.625,0.551)
#'dt=data.frame(t,i,f,fac)
#'fanan(dt)
fanan=function(dt)
{
  rw=nrow(dt)
    if (is.null(rw) || is.na(rw) || is.infinite(rw)) {
   stop(terminatepg)
    }
  terminatepg = "Terminating program ... \n"
  origmsg = "Here is the original message: "
  rw = tryCatch(rw)
  {
    as.integer(rw)
  }
  warning = function(warning_msg) {
    message("rw The  value is non-numeric, please try again")
    message(origmsg)
    message(paste(warning_msg, "\n"))
    return(NULL)
  }
  dt = tryCatch(dt)
    {
      as.data.frame(dt)
    }
    warning = function(warning_msg) {
      message(" The data frame is empty, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }

  sc=(2+dt[,1]-dt[,2]-dt[,3])/3
  ac=dt[,1]-dt[,3]
  ce=dt[,1]
  y1=sc
  y1=round(y1,2)
  y2=as.character(dt[,4])
  #------------------------
  ff=s_sort(y1,y2,ac,ce,rw)
  ff=s_sort(ac,y2,y1,ce,rw)
  ff=s_sort(ce,y2,ac,y1,rw)
  ff=s_sort(y1,y2,ac,ce,rw)
  y1=ff$y1
  y2=ff$y2
  ac=ff$ac
  ce=ff$ce
  #------------Processing------------
  for(i in 1: (rw-1))
  {
    if(y1[i]==y1[i+1])
    {
      if(ac[i]<ac[i+1])
      {
        k=ac[i]
        ac[i]=ac[i+1]
        ac[i+1]=k
        k=ce[i]
        ce[i]=ce[i+1]
        ce[i+1]=k
        k=y2[i]
        y2[i]=y2[i+1]
        y2[i+1]=k
      }
    }
  }
  #-----------------------
  y3=1
  for(i in 1:(rw-1))
  {
    if(y1[i]>y1[i+1])
    {
      k=y3[i]
      y3[i]=k
      y3[i+1]=i+1
    }
    if((y1[i]==y1[i+1]) & (ac[i]>ac[i+1]))
    {
      k=y3[i]
      y3[i]=k
      y3[i+1]=i+1
    }
    if((y1[i]==y1[i+1]) & (ac[i]<ac[i+1]))
    {
      k=y3[i]
      y3[i+1]=i+1
      y3[i]=k
    }
    if((y1[i]==y1[i+1]) & (ac[i]==ac[i+1]) & (ce[i]>ce[i+1] ))
    {
      k=y3[i]
      y3[i]=k
      y3[i+1]=i+1
    }
    if((y1[i]==y1[i+1]) & (ac[i]==ac[i+1]) & (ce[i]<ce[i+1] ))
    {
      k=y3[i]
      y3[i+1]=i+1
      y3[i]=k
    }
    if((y1[i]==y1[i+1]) & (ac[i]==ac[i+1]) & (ce[i]==ce[i+1] ))
    {
      y3[i]=(i+i+1)/2
      y3[i+1]= (i+i+1)/2
    }
  }
  sc=y1
  order=y3
  group=y2
  dt3=data.frame(sc,ac,ce,order,group)
  res=kruskal.test(order~group,dt3)
  message("the resulte of test","\n")
  message("----------------------------------","\n")
      message(paste0(capture.output(res),collapse="\n"))
      message("----------------------------------","\n")
  y2=as.integer(y2)
  nn=max(y2)
  rnk=matrix(nrow=nn)
  for(i in 1:nn)
  {
    rnk[i]=0
    for(j in 1:rw)
    {
      if(y2[j]==i)
      {rnk[i]=rnk[i]+order[j]}
    }
  }
  message("the essential data","\n")

    message("----------------------------------","\n")
    message(paste0(capture.output(dt),collapse="\n"))
    message("----------------------------------","\n")
    message("the data after operating","\n")
    message("----------------------------------","\n")
    message(paste0(capture.output(dt3),collapse="\n"))
    message("----------------------------------","\n")
    #----------------------- pair comparison----------------------
  for(j in 1:(nn-1))
    for(k in (j+1):nn)
    {
      {
        e=0
        s1=0
        s2=0
        for(i in 1:rw)
        {
          if(group[i]==j || group[i]==k)
          {e=e+1
          }
        }
        z1=matrix(nrow=e)
        z2=matrix(nrow=e)
        z3=matrix(nrow=e)
        r1=0
        r2=0
        s=0
        for(i in 1:rw)
        {
          if(group[i]==j || group[i]==k)
          {s=s+1
          z1[s]=sc[i]
          z2[s]=e-s+1
          z3[s]=group[i]}
        }
        for(i in 1:e)
        {
          if(z3[i]==j)
          {r1=r1+z2[i]
          s1=s1+1}
          if(z3[i]==k)
          {r2=r2+z2[i]
          s2=s2+1}
        }
        dt4=data.frame(z1,z2,z3)
        #print(dt4)
        #cat(r1,"    ",r2,"\n")
        un1=min(s1*s2+(s1*(s1+1)/2-r1), s1*s2+(s2*(s2+1)/2-r2))
        #print(un1)
        un2=s1*s2/2
        std=sqrt(s1*s2*(s1+s2+1)/12)
        zn=abs((un1-un2)/std)
        zn=round(zn,2)
        {message("z - calculated value =",zn,"\n")}

        if (zn>=1.96 )
        {
          if ((r1/s1)>=(r2/s2))
          {message(j ," better than ",k        ,"   *","\n")}
          else
          { message(k," better than ",j,"           *","\n")}
        }
        else
        { message("No difference between ",j," and ",k,"\n")}
        message("-------------------------------------------------------","\n")
      }}
}
