tryMe<-function(x){
	y1<-x^2
	y2<-x^3
	ret<-list(y1=y1, y2=y2)
	class(ret)<-"tryMe"
	return(ret)
}

print.tryMe<-function(x,...)
	print(x$y1)


