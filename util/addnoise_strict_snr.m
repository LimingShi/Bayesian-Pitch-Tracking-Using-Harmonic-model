function output_noise=addnoise_strict_snr(sig,input_noise,snr)
noise=input_noise;
noise_std_var=sqrt(10^(-snr/10)*(sig(:)'*sig(:))/(noise(:)'*noise(:)));
output_noise=noise_std_var*noise;
end
