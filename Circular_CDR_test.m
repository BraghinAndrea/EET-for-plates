function [hw,h_sigma,w_max,sigma_max] = EET_circolare_CDR_test(d,R2,R1,h_1,h_2,t,G,F)

%Restituisce i valori in [mm] degli spessori effettivi per deflessione e tensione massime 
% (risp. "hw" e "h_sigma"), la deflessione massima al centro della piastra in [mm] ("w_max")
% e la tensione massima in [Mpa] ("sigma_max").

%Le grandezze di input richieste sono il diametro della piastra circolare
%("d"), la distanza radiale dei punti di appoggio della piastra ("R2"), la 
% Distanza radiale dei punti di applicazione del carico ("R1")
% gli spessori dei due strati di vetro in [mm] ("h_1" e "h_2"), lo spessore dell'intercalare in [mm] ("t"), 
% il modulo a taglio dell'intercalare in [Mpa] ("G") ed il carico in [N/mm] ("F");

% NB: R3>R2>R1

R3=d/2;

nu_glass=0.22; %coefficiente di Poisson del vetro 
E_glass=70*10^3; %Modulo elastico del vetro [MPa]

psi=(16*pi*(nu_glass + 1)*(301*nu_glass - (1800)*log(2) + (900)*log(3) - 1800*nu_glass*log(2) + 900*nu_glass*log(3) + 399))/(R3^2*pi*1160*log(2) - 462*nu_glass - (972)*log(3) + (392)*log(5) + 2320*nu_glass*log(2) - 1944*nu_glass*log(3) + 784*nu_glass*log(5) + 1160*nu_glass^2*log(2) - 972*nu_glass^2*log(3) + 392*nu_glass^2*log(5) - 329*nu_glass^2 + 259);

H=t+((h_1+h_2)/2);

hs_1=H*h_1/(h_1+h_2);
hs_2=H*h_2/(h_1+h_2);

Is=h_1*hs_2^2+h_2*hs_1^2;

D_1=E_glass*h_1^3/(12*(1-nu_glass^2));
D_2=E_glass*h_2^3/(12*(1-nu_glass^2));
D_tot=D_1+D_2+(12*D_1*D_2*H^2/((D_1*h_2^2)+(D_2*h_1^2)));

eta=1./(1.+((t./G).*((D_1+D_2)/D_tot)*(12*D_1*D_2*psi)/((D_1*h_2^2)+(D_2*h_1^2))));
D_R=((eta./D_tot)+((1-eta)./(D_1+D_2))).^(-1); %rigidezza effettiva equivalente

hw=(12.*D_R.*(1-nu_glass^2)./E_glass).^(1/3);
h_sigma=1./((2.*eta*hs_1/(h_1^3+h_2^3+12*Is))+(h_2./hw.^3)).^0.5;

w_max=(F*((18*R3^4*log((3*R3)/5))/25 - (18*R3^4*log((4*R3)/5))/25 + (238*R3^4*nu_glass)/625 + (462*R3^4)/625 + (18*R3^4*nu_glass*log((3*R3)/5))/25 - (18*R3^4*nu_glass*log((4*R3)/5))/25))./(16.*D_R.*pi.*R3^2.*(nu_glass + 1));

sigma_max=(3*F./(4*pi.*h_sigma.^2)).*(((1-nu_glass)*(1-(R1^2/R2^2))*(R2^2/R3^2))+(2*(1+nu_glass)*log(R2/R1)));

end