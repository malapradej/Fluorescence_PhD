function Fout = meanleaf(canopy,F,choice,Ps)

nl    = canopy.nlayers;
nli   = canopy.nlincl;
nlazi = canopy.nlazi;
lidf  = canopy.lidf;

% author: Dr. ir. Christiaan van der Tol (tol@itc.nl)
% date:     7   December 2007
% update:   11  February 2008 made modular (Joris Timmermans)

% update:   25 Feb 2013 Wout Verhoef : Propose name change, remove globals
%                                      and use canopy-structure for input
%
% function [F2,F3] = F1tot(F,choice,Ps)
% calculates the layer average and the canopy average of leaf properties
% per layer, per leaf angle and per leaf azimuth (36)
%
% Input:
%   F       input matrix (3D)   [nli, nlazi,nl]
%   choice  integration method  'angles'            : integration over leaf angles
%                               'angles_and_layers' : integration over leaf layers and leaf angles
%   Ps      fraction sunlit per layer [nl]
%
% Output:
%   Fout    in case of choice = 'angles': [nl]
%           in case of choice = 'angles_and_layers': [1]
%
% J.Malaprade: amended the function to include integration over angles but 
% conserving the spectral information. Use 'angles_spectral' as choice.
% Also included integration over angles and layers but conserving spectral
% information. Use 'angles_layers_spectral' as choice.

if strcmp(choice, 'angles_spectral') || strcmp(choice, 'angles_layers_spectral')
    nwlPARFs = length(canopy.IparFs);
    Fout = zeros(nli, nlazi, nl, nwlPARFs);
else
    Fout = zeros(nli, nlazi,nl);
end

switch choice
%% integration over leaf angles
    case 'angles'
        
        for j = 1:nli                                   
             Fout(j,:,:)    = F(j,:,:)*lidf(j);         % [nli, nlazi,nl]
        end
        Fout                = sum(sum(Fout))/nlazi;     % [1,1,nl]
        Fout                = permute(Fout,[3 1 2]);    % [nl]
%% integration over leaf angles keeping spectral information
    case 'angles_spectral'
        
        for j = 1:nli                                   
             Fout(j,:,:,:)    = F(j,:,:,:)*lidf(j);         % [nli,nlazi,nl,nwlPARFs]
        end
        Fout                = sum(sum(Fout))/nlazi;     % [1,1,nl,nwlPARFs]
        Fout                = permute(Fout,[3 1 2 4]);    % [nl,nwlPARFs]
%% integration over leaf angles and layers keeping spectral information
    case 'angles_layers_spectral'
        
        for j = 1:nli
            Fout(j,:,:,:)   = F(j,:,:,:)*lidf(j);           % [nli,nlazi,nl,nwlPARFs]
        end
        
        for j = 1:nl
            Fout(:,:,j,:)   = Fout(:,:,j,:)*Ps(j);
        end
        Fout                = sum(sum(sum(Fout)))/nlazi/nl; % [nwlPARFs]
%% integration over layers only         
    case 'layers'
        %not implemented
%% integration over both leaf angles and layers       
    case 'angles_and_layers'
        for j = 1:nli
             Fout(j,:,:)      = F(j,:,:)*lidf(j);
        end
        
        for j = 1:nl
             Fout(:,:,j)    = Fout(:,:,j)*Ps(j);
        end
        Fout                = sum(sum(sum(Fout)))/nlazi/nl;
end