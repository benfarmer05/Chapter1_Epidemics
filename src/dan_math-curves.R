x = seq(0,1,.001) 
a = 10^(seq(-2,1,.2)) 
y = x
plot(y~x,type = "n") 

for (i in 1:length(a)) {   
  y = x^a[i]   
  lines(y~x) 
}
