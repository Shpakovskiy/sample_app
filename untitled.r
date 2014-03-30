

N=140	#numbers of particles
T=80	#time of calculation
zero <- array(rep(0,N*2*T), dim=c(N,2,T))

a=100	#side of squer
u=10	#velocity of particles
tau=0.1*a/u
sigma=tau*u

init_coordinate <-function(coordinate, a)	#initial coordinate
  {
  N <- dim(coordinate)[1]
  for (n in 1:N)
    {
    coordinate[n,1,1] <- a/2
    coordinate[n,2,1] <- 3*a/2
    }
  return (coordinate)
  }
  
init_velocity <-function(velocity, u)	#initial velocity
  {
  N <- dim(velocity)[1]
  for (n in 1:N)
    {
    velocity[n,1,1] <- u*sin(n/N*2*pi)
    velocity[n,2,1] <- u*cos(n/N*2*pi)
    }
  return (velocity)
  }

velocity <- init_velocity(zero, u)
coordinate <- init_coordinate(zero, a)

inter_x <-function(x,y,a,sigma)	#interaction V_x with walls
  {
  v=1
  if ( (y>a & (x<0+sigma | x>7*a-sigma)) | (y<a & (x<a+sigma | x>8*a-sigma)) )
      {
      v <--v
      }
  return (v)
  }
  
inter_y <-function(x,y,a,sigma)	#interaction V_y with walls
  {
  v=1
  if ( y>2*a-sigma | y<0+sigma | (y<a+sigma & x<a+sigma) | (y>a-sigma & x>6*a-sigma ) )
      {
      v <--v
      }
  return (v)
  }
  

lagrange <-function(coordinate,velocity,T,a,u,tau)	#calculation of distribution
  {
  T=dim(coordinate)[3]
  N=dim(coordinate)[1]
  k=2	#numbers of squers
  p=8	#numbers of squers
  sq <- array(rep(0,k*p*T), dim=c(k,p,T))
  stop <- array(rep(0,k*p), dim=c(k,p))
  relaks <- array(rep(0,k*p), dim=c(k,p))
  dispersion <- array(rep(0,k*p*T), dim=c(k,p,T))
  for (t in 2:T)
    {
    
    for (n in 1:N)
      {
      
      x <- coordinate[n,1,t-1]
      y <- coordinate[n,2,t-1]
      u_x <- velocity[n,1,t-1]*inter_x(x,y,a,sigma)
      u_y <- velocity[n,2,t-1]*inter_y(x,y,a,sigma)
      coordinate[n,1,t] <- x+u_x*tau
      coordinate[n,2,t] <- y+u_y*tau
      velocity[n,1,t] <-u_x
      velocity[n,2,t] <-u_y
      
      for (i in 1:k)
	{
	for (j in 1:p)
	  {
	  if ((x > a*(j-1)) & (x < a*j) & (y > a*(i-1)) & (y < a*i))
	    {
	    sq[i,j,t] <- sq[i,j,t]+1
	    }
	  dispersion[i,j,t] <- round(abs(sq[i,j,t]-N/(p*k-2))/N, digits=4)
	  if ((dispersion[i,j,t]<0.01) & (stop[i,j] != 1))
	    {
	    relaks[i,j] <- t
	    stop[i,j] <- 1
	    }
	  }
	}
      }
    
    }
  return (dispersion)
  }
coordinate <- lagrange(coordinate,velocity,T,a,u,tau)
print (coordinate)
