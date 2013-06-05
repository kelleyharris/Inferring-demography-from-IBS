
######
#You need to set theta and rho such that theta = 40000*mu, where mu is the mutation rate per site per generation, and rho = 40000*r, where r is the recombination rate per site per generation
######
global theta
theta=0.001 
global rho
rho=0.0004
#######

import math
from math import log
from math import exp
from math import pi

from mpmath import gamma
from mpmath import gammainc

import cmath

from mpmath import hyp1f1 as hyper

def prob_L_from_mut(L,ta,ts,phi,N, theta0=theta):
    if L<200000 and L>0:
        if phi>0.01:
            b=0.5+0.25*L*rho*N+0.5*L*theta*N
            c=(0.25+0.5*exp(ta/N)-0.5*exp(2*ta/N))*L*rho*N
            prob1=phi*0.5*exp(0.5*L*((-1.5+exp(-ta/N))*N*rho-2*ta*theta)+c)*hyper(1,b+1,-c)/b
        #        prob=phi*0.5*exp(0.5*L*((-1.5+exp(-ta/N))*N*rho-2*ta*theta))*(-c)**(-b)*(gamma(1+b)-b*gammainc(b,0,-c))/b
            b+=0.5*theta0*N
            c=(0.25+exp(ta/N)*(0.5-0.5*exp(ta/N)))*L*rho*N
            prob2=-phi*0.5*exp(-0.75*L*rho*N+0.5*exp(-ta/N)*L*rho*N-ta*theta0-L*ta*theta+c)*hyper(1,b+1,-c)/b
        else:
            prob1, prob2 = 0,0
        if phi<0.99:
            b=0.5+0.25*L*rho*N+0.5*L*theta*N
            c=(0.25+0.5*exp(ts/N)-0.5*exp(2*ts/N))*L*rho*N
        #        prob+=(1-phi)*0.5*exp(0.5*L*((-1.5+exp(-ts/N))*N*rho-2*ts*theta))*(-c)**(-b)*(gamma(1+b)-b*gammainc(b,0,-c))/b
            prob3=(1-phi)*0.5*exp(0.5*L*((-1.5+exp(-ts/N))*N*rho-2*ts*theta)+c)*hyper(1,b+1,-c)/b
            b+=0.5*theta0*N
            c=(0.25+exp(ts/N)*(0.5-0.5*exp(ts/N)))*L*rho*N
            prob4=-(1-phi)*0.5*exp(-0.75*L*rho*N+0.5*exp(-ts/N)*L*rho*N-ts*theta0-L*ts*theta+c)*hyper(1,b+1,-c)/b
        else:
            prob3, prob4=0,0
        return abs(math.fsum((prob1,prob2,prob3,prob4)))

    else:
        prob=phi*(math.exp(-ta*L*(rho+theta))/(1+L*N*(rho+theta))-math.exp(-ta*(L*(rho+theta)+theta))/(1+L*(rho+theta)*N+theta*N))*(1+L*((-exp(-rho*ta)+1)/(1-rho*N)-rho*N*(1-exp(-ta))/(1-rho*N)))
        prob+=(1-phi)*(math.exp(-ts*L*(rho+theta))/(1+L*(rho+theta)*N)-math.exp(-ts*(L*(rho+theta)+theta))/(1+L*(rho+theta)*N+theta*N))*(1+L*((-exp(-rho*ts)+1)/(1-rho*N)-rho*N*(1-exp(-ts))/(1-rho*N)))
        return prob

def recent_shallow_recent(L,phi,ta,ts,N):
    prob=math.log((1+(L-1)*(rho+theta)*N-rho*N)*(1+(L-1)*N*(rho+theta))/((1+theta*N)*(1+rho*N+theta*N)))/(2+L*(rho+theta)*N-rho*N)
    prob-=math.exp(-ta*theta)*math.log((1+(L-1)*N*(rho+theta))*(1+(L-1)*N*(rho+theta))/((1+N*theta)*(1+N*(rho+2*theta))))/(2+(L*(rho+theta)+theta-rho)*N)
    prob*=phi*phi*math.exp(-ta*L*(rho+theta))*(1-math.exp(-ta*rho))/((rho+theta)*N)
    return prob

def recent_deep_recent(L,phi,ta,ts,N):
    prob=math.log((1+(L-1)*(rho+theta)*N)/(1+rho*N+theta*N))*math.exp(-ta*L*(rho+theta))*(1.0/(1+L*N*(rho+theta))-1.0/(2+L*N*(rho+theta)-rho*N))
    prob-=math.log((1+L*(rho+theta)*N-rho*N)/(1+rho*N+2*theta*N))*math.exp(-ta*(L*(rho+theta)+theta))*(1.0/(1+L*(rho+theta)*N+theta*N)-1.0/(2+L*(rho+theta)*N+N*theta-N*rho))
    prob+=math.log((1+(L-1)*(rho+theta)*N-rho*N)/(1+theta*N))*math.exp(-ta*L*(rho+theta))*(1.0/(1+L*(rho+theta)*N)-1.0/(2+L*(rho+theta)*N-rho*N)-math.exp(-ta*theta)/(1+L*(rho+theta)*N+theta*N)+math.exp(-ta*theta)/(2+L*(rho+theta)*N+theta*N-rho*N))
    return prob*phi*rho/((1-rho*N)*(rho+theta))

def ancient_deep_ancient(L,phi,ta,ts,N):
    prob=math.log((1+(L-1)*(rho+theta)*N)/(1+rho*N+theta*N))*math.exp(-ts*L*(rho+theta))*(1.0/(1+L*N*(rho+theta))-1.0/(2+L*N*(rho+theta)-rho*N))
    prob-=math.log((1+L*(rho+theta)*N-rho*N)/(1+rho*N+2*theta*N))*math.exp(-ts*(L*(rho+theta)+theta))*(1.0/(1+L*(rho+theta)*N+theta*N)-1.0/(2+L*(rho+theta)*N+N*theta-N*rho))
    prob+=math.log((1+(L-1)*(rho+theta)*N-rho*N)/(1+theta*N))*math.exp(-ts*L*(rho+theta))*(1.0/(1+L*(rho+theta)*N)-1.0/(2+L*(rho+theta)*N-rho*N)-math.exp(-ts*theta)/(1+L*(rho+theta)*N+theta*N)+math.exp(-ts*theta)/(2+L*(rho+theta)*N+theta*N-rho*N))
    return prob*(1-phi)*rho/((1-rho*N)*(rho+theta))
 

def ancient_shallow_ancient(L,phi,ta,ts,N):
    prob=math.log((1+(L-1)*(rho+theta)*N-rho*N)*(1+(L-1)*N*(rho+theta))/((1+theta*N)*(1+rho*N+theta*N)))/(2+L*(rho+theta)*N-rho*N)
    prob-=math.exp(-ts*theta)*math.log((1+(L-1)*N*(rho+theta))*(1+(L-1)*N*(rho+theta))/((1+N*theta)*(1+N*(rho+2*theta))))/(2+(L*(rho+theta)+theta-rho)*N)
    prob*=(1-phi)*(1-phi)*math.exp(-ts*L*(rho+theta))*(1-math.exp(-ts*rho))/((rho+theta)*N)
    return prob

def ancient_medium_ancient(L,phi,ta,ts,N):
    prob=math.log((1+N*(L-1)*(rho+theta)-N*rho)*(1+N*(L-1)*(rho+theta))/((1+N*theta)*(1+N*rho+N*theta)))/(2+L*(rho+theta)*N-rho*N)
    prob-=math.exp(-ts*theta)*math.log((1+(L-1)*(rho+theta)*N)*(1+(L-1)*(rho+theta)*N)/((1+theta*N)*(1+rho*N+2*theta*N)))/(2+L*(rho+theta)*N+theta*N-rho*N)
    prob*=(1-phi)*math.exp(-ts*L*(rho+theta))*(math.exp(-ta*rho)-math.exp(-ts*rho))/(N*(rho+theta))
    return prob

def recent_shallow_ancient(L,phi,ta,ts,N):
    prob1=math.exp(-(ts-ta)*theta)*math.log(1+1.0/((ts-ta)*(1+theta)))-0.5*math.exp(-(ts-ta)*((L-1)*(rho+theta)-rho))*math.log(1+2.0/((ts-ta)*(1+(L-1)*(rho+theta)-rho)))
    prob1+=math.exp(-(ts-ta)*theta)*math.log(1+1.0/((ts-ta)*(1+(L-1)*(rho+theta))))-0.5*math.exp(-(ts-ta)*((L-1)*(rho+theta)-rho))*math.log(1+2.0/((ts-ta)*(1+rho+theta)))
    prob1*=math.exp(-ta*(L*(rho+theta)-rho))/(2+L*(rho+theta)-rho)
    prob2=math.exp(-(ts-ta)*theta)*math.log(1+1.0/((ts-ta)*(1+theta)))-0.5*math.exp(-(ts-ta)*((L-1)*(rho+theta)-rho))*math.log(1+2.0/((ts-ta)*(1+(L-1)*(rho+theta)-rho)))
    prob2+=math.exp(-(ts-ta)*theta)*math.log(1+1.0/((ts-ta)*(1+L*(rho+theta)-rho)))-0.5*math.exp(-(ts-ta)*((L-1)*(rho+theta)-rho))*math.log(1+2.0/((ts-ta)*(1+2*theta+rho)))
    prob2*=math.exp(-ta*(L*(rho+theta)+theta-rho))/(2+L*(rho+theta)+theta-rho)
    return (prob1-prob2)*phi*(1-phi)*(1-math.exp(-ta*rho))/(rho+theta)


def ancient_shallow_recent(L,phi,ta,ts,N):
    prob1=math.exp(-(ts-ta)*(rho+theta))*math.log(1+1.0/((ts-ta)*(1+rho+theta)))-0.5*math.exp(-(ts-ta)*(L-1)*(rho+theta))*math.log(1+2.0/((ts-ta)*(1+(L-1)*(rho+theta))))
    prob1+=math.exp(-(ts-ta)*(rho+theta))*math.log(1+1.0/((ts-ta)*(1+(L-1)*(rho+theta)-rho)))-0.5*math.exp(-(ts-ta)*((L-1)*(rho+theta)))*math.log(1+2.0/((ts-ta)*(1+theta)))
    prob1*=math.exp(-ta*(L*(rho+theta)-rho))/(2+L*(rho+theta)-rho)
    prob2=math.exp(-(ts-ta)*(2*theta+rho))*math.log(1+1.0/((ts-ta)*(1+2*theta+rho)))-0.5*math.exp(-(ts-ta)*(L*(rho+theta)-rho))*math.log(1+2.0/((ts-ta)*(1+L*(rho+theta)-rho)))
    prob2+=math.exp(-(ts-ta)*(2*theta+rho))*math.log(1+1.0/((ts-ta)*(1+(L-1)*(rho+theta)-rho)))-0.5*math.exp(-(ts-ta)*(L*(rho+theta)-rho))*math.log(1+2.0/((ts-ta)*(1+theta)))
    prob2*=math.exp(-ta*(L*(rho+theta)+theta-rho))/(2+L*(rho+theta)+theta-rho)
    return (prob1-prob2)*phi*(1-phi)*(1-math.exp(-ta*rho))/(rho+theta)

def prob_L_from_mut_approx(L,ta,ts,phi,N):
    prob=phi*(math.exp(-ta*L*(rho+theta))/(1+L*(rho+theta))-math.exp(-ta*(L+1)*(rho+theta))/(1+(L+1)*(rho+theta)))
    prob+=(1-phi)*(math.exp(-ts*L*(rho+theta))/(1+L*(rho+theta))-math.exp(-ts*(L+1)*(rho+theta))/(1+(L+1)*(rho+theta)))
    prob+=s1_minus_s2_approx(L,N)*(math.exp(ta*rho)-1)*(phi*phi*math.exp(-ta*L*(rho+theta))+(1-phi)*(1-phi)*math.exp(-ts*L*(rho+theta)))
    prob+=s2_approx(L,ta,ts,phi,N)*(math.exp(ta*rho)-1)*(phi*phi*math.exp(-ta*L*(rho+theta))*ta*(rho+theta)+(1-phi)*(1-phi)*math.exp(-ts*L*(rho+theta))*ts*(rho+theta))
#    prob+=s1(L,ta,ts,phi,N)*((phi*phi*math.exp(-ta*L*N*(rho+theta))+(1-phi)*(1-phi)*math.exp(-ts*N*L*(rho+theta)))*(math.exp(ta*rho*N)-1)+(1-phi)*math.exp(-ts*L*N*(rho+theta))*(math.exp(rho*N*(ts-ta))-1))
#    prob-=s2(L,ta,ts,phi,N)*((phi*phi*math.exp(-ta*(L+1)*N*(rho+theta))+(1-phi)*(1-phi)*math.exp(-ts*N*(L+1)*(rho+theta)))*(math.exp(ta*rho*N)-1)+(1-phi)*math.exp(-ts*(L+1)*N*(rho+theta))*(math.exp(rho*N*(ts-ta))-1))
    prob+=N*rho*(phi*math.exp(-ta*L*N*(rho+theta))+(1-phi)*math.exp(-ts*L*N*(rho+theta)))*s7_8_approx(L,ta,ts,phi,N)
    prob-=N*rho*(phi*math.exp(-ta*(L+1)*N*(rho+theta))+(1-phi)*math.exp(-ts*N*(L+1)*(rho+theta)))*s9_10_approx(L,ta,ts,phi,N)
    return prob

## This defines a cubic that allows for a C^1 elementary approximation to the polylog

s1=log(4)/4
y0=-pi**2/6-log(4)**2/2

s2=log(2)
y1=-pi**2/12-log(2)

x0=-4
x1=-2

global a_poly
a_poly=-(-s1*x0-s2*x0+s1*x1+s2*x1+2*y0-2*y1)/(x0-x1)**3
global b_poly
b_poly=-(s1*x0**2+2*s2*x0**2+s1*x0*x1-s2*x0*x1-2*s1*x1**2-s2*x1**2-3*x0*y0-3*x1*y0+3*x0*y1+3*x1*y1)/(x0-x1)**3
global c_poly
c_poly=-(-s2*x0+s1*x1)/(x0-x1)-3*x0*x1*(-s1*x0-s2*x0+s1*x1+s2*x1+2*y0-2*y1)/(x0-x1)**3
global d_poly
d_poly=-(s2*x0**3*x1+s1*x0**2*x1**2-s2*x0**2*x1**2-s1*x0*x1**3-3*x0*x1**2*y0+x1**3*y0-x0**3*y1+3*x0**2*x1*y1)/(x0-x1)**3

def Li(z):
    if z>-2:
        return -pi**2/12+(z+1)*log(2)
    elif z>-4:
        return a_poly*z**3+b_poly*z**2+c_poly*z+d_poly
    else:
        return -pi**2/6-log(-z)**2/2

def prob_L_2_recombs(L,ts,N):
    t=theta*N
    r=rho*N
    prob1=0
    prob=0

    prob1+=(log((r+(L-2)*(r+t))/(2*r+t))*log(1+(L-1)*(r+t))-Li((r+(L-2)*(r+t))/(1+(L-1)*(r+t)))+Li((2*r+t)/(1+(L-1)*(r+t))))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t)))
    prob1-=log(1+t)*log((r+(L-2)*(r+t))/(2*r+t))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t))) ## integrate 1/(rho+L1(rho+theta))
    prob1-=(log(2+t+(L-2)*(r+t))*log((1+(L-2)*(r+t))/(1+r+t))-Li((1+(L-2)*(r+t))/(2+t+(L-2)*(r+t)))+Li((1+r+t)/(2+t+(L-2)*(r+t))))/((1-r)*(3-2*r+L*(r+t))*(2-r+L*(r+t)))
    prob1+=log(1+t)*log((1+(L-2)*(r+t))/(1+r+t))/((1-r)*(3-2*r+L*(r+t))*(2-r+L*(r+t))) ## integrate 1/(1+L1(rho+theta))
    prob1+=(-log(1+t)*log((1+r+2*t)/(r+t))+log(1+t+(L-3)*(r+t))*log((1+t+(L-2)*(r+t))/(r+t))-Li(-(1+t)/(r+t))+Li(-(1+t+(L-3)*(r+t))/(r+t)))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t)))
    prob1-=log(1+t)*log((1-r+(L-1)*(r+t))/(1+2*r+t))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t))) ## integrate 1/(1-rho+(L-L1)(rho+theta))
    prob1-=(Li(-(1+t+(L-3)*(r+t))/(1+t))+log(1+t+(L-3)*(r+t))*log((2-2*r+(L-1)*(r+t))/(1+t))-Li(-1)-log(1+t)*log(2))/((1-r)*(3-2*r+L*(r+t))*(2-r+L*(r+t)))
    prob1+=log(1+t)*log((2-2*r+(L-1)*(r+t))/(2+2*t))/((1-r)*(3-2*r+L*(r+t))*(2-r+L*(r+t))) ##integrate 1/(2-2rho+(L-L1)(rho+theta))

    prob2=-(log(2+t+(L-2)*(r+t))*log((1+(L-2)*(r+t))/(1+r+t))-Li((1+(L-2)*(r+t))/(2+t+(L-2)*(r+t)))+Li((1+r+t)/(2+t+(L-2)*(r+t))))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t)))
    prob2+=log(1+t)*log((1+(L-2)*(r+t))/(1+r+t))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t))) ## integrate 1/(1+L1(rho+theta))
    prob2+=(log(1+(L-1)*(r+t))*log((r+(L-2)*(r+t))/(2*r+t))+Li((2*r+t)/(1+(L-1)*(r+t)))-Li((r+(L-2)*(r+t))/(1+(L-1)*(r+t))))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t)))
    prob2-=log(1+t)*log((r+(L-2)*(r+t))/(2*r+t))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t))) ## integrate 1/(rho+L1(rho+theta))

    prob3=log(1+(L-1)*(r+t))*log((r+(L-2)*(r+t))/(2*r+t))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t)))
    prob3-=(log(1+t)*log((r+(L-2)*(r+t))/(2*r+t))-Li(-(r+(L-2)*(r+t))/(1+t))+Li(-(2*r+t)/(1+t)))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t))) ## integrate 1/(rho+L1(rho+theta))
    prob3-=log(1+(L-1)*(r+t))*log((1+(L-2)*(r+t))/(1+r+t))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t)))
    prob3+=(log(r+t)*log((1+(L-2)*(r+t))/(1+r+t))-Li(-(1+(L-2)*(r+t))/(r+t))+Li(-(1+r+t)/(r+t)))/((1-r)*(1+L*(r+t))*(2-r+L*(r+t)))  ## integrate 1/(1+L1(rho+theta))

    prob=r**2/(r+t)**2*exp(-ts*L*(rho+theta))*(2*prob1-prob2-prob3)    

    prob4=0
    prob4=(log((L-1)/2)*log(1+t+(L-1)*(r+t))-Li((L-1)*(r+t)/(1+t+(L-1)*(r+t)))+Li(2*(r+t)/(1+t+(L-1)*(r+t))))/((1-r)*(1+t+L*(r+t))*(2-r+t+L*(r+t)))
    prob4-=log(1+t)*log((L-1)/2)/((1-r)*(1+t+L*(r+t))*(2-r+t+L*(r+t))) ## integrate 1/(1+L1)
    prob4-=(log(2+2*t+(L-2)*(r+t))*log((1+t+(L-2)*(r+t))/(1+2*t+r))-Li((1+t+(L-2)*(r+t))/(2+2*t+(L-2)*(r+t)))+Li((1+2*t+r)/(2+2*t+(L-2)*(r+t))))/((1-r)*(3-2*r+t+L*(r+t))*(2-r+t+L*(r+t)))
    prob4+=log(1+t)*log((1+t+(L-2)*(r+t))/(1+2*t+r))/((1-r)*(3-2*r+t+L*(r+t))*(2-r+t+L*(r+t))) ## integrate 1/(1+theta+L1(rho+theta))
    prob4+=(-log(1+t)*log((1+r+2*t)/(r+t))+log(1+t+(L-3)*(r+t))*log((1+t+(L-2)*(r+t))/(r+t))-Li(-(1+t)/(r+t))+Li(-(1+t+(L-3)*(r+t))/(r+t)))/((1-r)*(1+t+L*(r+t))*(2-r+t+L*(r+t)))
    prob4-=log(1+t)*log((1-r+(L-1)*(r+t))/(1+2*r+t))/((1-r)*(1+t+L*(r+t))*(2-r+t+L*(r+t))) ## integrate 1/(1-rho+(L-L1)(rho+theta))
    prob4-=(Li(-(1+t+(L-3)*(r+t))/(1+t))+log(1+t+(L-3)*(r+t))*log((2-2*r+(L-1)*(r+t))/(1+t))-Li(-1)-log(1+t)*log(2))/((1-r)*(3-2*r+t+L*(r+t))*(2-r+t+L*(r+t)))
    prob4+=log(1+t)*log((2-2*r+(L-1)*(r+t))/(2+2*t))/((1-r)*(3-2*r+t+L*(r+t))*(2-r+t+L*(r+t))) ## integrate 1/(2-2*rho+(L-L1)*(rho+theta))
    
    prob5=-(log(2+2*t+(L-2)*(r+t))*log((1+t+(L-2)*(r+t))/(1+2*t+r))-Li((1+t+(L-2)*(r+t))/(2+2*t+(L-2)*(r+t)))+Li((1+2*t+r)/(2+2*t+(L-2)*(r+t))))/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t)))
    prob5+=log(1+t)*log((1+t+(L-2)*(r+t))/(1+2*t+r))/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t)))
    prob5+=(log(1+t+(L-1)*(r+t))*log((L-1)/2)+Li(2*(r+t)/(1+t+(L-1)*(r+t)))-Li((L-1)*(r+t)/(1+t+(L-1)*(r+t))))/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t))) ##changed one L-2 to L-1 and it made things slightly worse
    prob5-=log(1+t)*log((L-1)/2)/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t)))
    
    prob6=log(1+t+(L-1)*(r+t))*log((L-1)/2)/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t)))
    prob6-=(log(1+t)*log((L-1)/2)-Li(-(L-1)*(r+t)/(1+t))+Li(-2*(r+t)/(1+t)))/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t)))
    prob6-=log(1+t+(L-1)*(r+t))*log((1+t+(L-2)*(r+t))/(1+2*t+r))/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t)))
    prob6+=(log(r+t)*log((1+t+(L-2)*(r+t))/(1+2*t+r))-Li(-(1+t+(L-2)*(r+t))/(r+t))+Li(-(1+r+2*t)/(r+t)))/((1-r)*(2-r+t+L*(r+t))*(1+t+L*(r+t)))
    
    prob-=r**2/(r+t)**2*exp(-ts*(L*(rho+theta)+theta))*(2*prob4-prob5-prob6)


    # three recent recombinations:(bug-checked)

    prob13=1.0/((r+t)*(3-2*r+L*(r+t)))*(log((2+(-3+L)*r+(-2+L)*t)/(2+t))*log((1+(-2+L)*r+(-1+L)*t)*(1+(-1+L)*(r+t)))+log((1+(-2+L)*r+(-1+L)*t)/(1+r+2*t))*log((1+(-2+L)*r+(-1+L)*t)*(1+(-1+L)*(r+t))))
    prob13-=log((1+t)*(1+r+t))*(-log(2+t)-log(1+r+2*t)+log(2+(-3+L)*r+(-2+L)*t)+log(1+(-2+L)*r+(-1+L)*t))/((r+t)*(3-2*r+L*(r+t)))
    prob13-=exp(-ts*theta)*2.0/((r+t)*(3-2*r+t+L*(r+t)))*(log(1+(-3+L)*(r+t))*log((2+(-3+L)*r+(-1+L)*t)/(1+2*t))-Li(-1.0/(1+2*t))+Li(-(1+(-3+L)*(r+t))/(1+2*t))+(-log(1+r+2*t)+log(1+(-2+L)*r+(-1+L)*t))*log(1+(-1+L)*(r+t)))
    prob13+=exp(-ts*theta)*log(1+t)*log(1+r+2*t)*(-log(2*(1+t))-log(1+r+2*t)+log(2+(-3+L)*r+(-1+L)*t)+log(1+(-2+L)*r+(-1+L)*t))/((r+t)*(3+(-2+L)*r+t+L*t))
    
    prob+=prob13*(1-exp(-ts*rho))**2*exp(-ts*(L*(rho+theta)-2*rho))


    return prob
    

def prob_L_from_mut_precise_varmu(L,ta,ts,phi,N, theta0=theta):
    prob=0.999*prob_L_from_mut(L,ta,ts,phi,N)
    prob+=0.001*prob_L_from_mut(L,ta,ts,phi,N, theta0)    
    if phi>0.01:
        prob+=recent_shallow_recent(L,phi,ta,ts,N)
        prob+=recent_deep_recent(L,phi,ta,ts,N)
    if phi<0.99:
        prob+=ancient_shallow_ancient(L,phi,ta,ts,N)
#    prob+=recent_shallow_ancient(L,phi,ta,ts,N)
#    prob+=2*ancient_shallow_recent(L,phi,ta,ts,N)
        prob+=ancient_deep_ancient(L,phi,ta,ts,N)
#        prob+=ancient_medium_ancient(L,phi,ta,ts,N)
    if L>3:
        if phi>0.01:
            prob+=phi*prob_L_2_recombs(L,ta,N)
        if phi<0.99:
            prob+=(1-phi)*prob_L_2_recombs(L,ts,N)
    return prob
    
