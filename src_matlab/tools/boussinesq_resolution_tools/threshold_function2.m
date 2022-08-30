function OUT=threshold_function2(x)
% sigmoid function centered on 1 so it is concave
% %     OUT=1-1./(1+exp(-100*(x-0.99)));
%     OUT=1-2./(1+exp(-1000*(x-1)));
    OUT=1-exp(-10*(1-x));
%     OUT=double(x>1)+(x<=1).*(1-exp(-1000*(1-x)));
% % % % %     OUT=1-1.01./(1+exp(-((x-1))/1e-2+log(0.01)));%2./(1+exp(-(x-1)/1e-3));
% % % % % % % % % %     K=1.0001;
% % % % % % % % % %     nu=0.01;%100;%5;
% % % % % % % % % %     B=100;%1000;
% % % % % % % % % %     OUT=1-K./(1+(K^nu-1)*exp(-B*(x-1))).^(1/nu);
%     OUT=1-exp(-100*(1-x));
%     OUT=1*ones(size(x));
%     OUT=1-1./(1+exp(-1000*(x-0.99)));   
% OUT=x<1;
%     OUT=(1-1./(1+exp(-600*(x-0.97))));
%     OUT=(1-1./(1+exp(-70*(x-0.8))));
end