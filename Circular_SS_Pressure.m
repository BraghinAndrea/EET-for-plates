function [hw,h_sigma,w_max,sigma_max] = EET_circolare_appoggiata_carico_distribuito(d,h_1,h_2,t,G,q)

%Restituisce i valori in [mm] degli spessori effettivi per deflessione e tensione massime 
% (risp. "hw" e "h_sigma"), la deflessione massima al centro della piastra in [mm] ("w_max")
% e la tensione massima in [Mpa] ("sigma_max").

%Le grandezze di input richieste sono il diametro della piastra circolare ("d"),
% gli spessori dei due strati di vetro in [mm] ("h_1" e "h_2"), lo spessore dell'intercalare in [mm] ("t"), 
% il modulo a taglio dell'intercalare in [Mpa] ("G") ed il carico in [N/mm^2] ("q");

R=d/2; %Raggio della piastra circolare
nu_glass=0.22; %coefficiente di Poisson del vetro 
E_glass=70*10^3; %Modulo elastico del vetro [MPa]

psi=(16*(7+nu_glass)*(1+nu_glass))/(R^2*(nu_glass^2+10*nu_glass+33));

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

w_max=(5+nu_glass)*q*R^4./(64.*D_R.*(1+nu_glass));
sigma_max=3*q*(R^2)*(3+nu_glass)./(8.*h_sigma.^2);

end