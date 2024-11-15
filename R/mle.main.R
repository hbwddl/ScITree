mle.main<-function(){
  cpp_out<-mle_cpp()
  cat("\nRequested likelihood generated successfully.\n ")
  return(cpp_out)
}