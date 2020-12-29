clear all; clc;
clf;

% === plot the energy change with respect to time

[time, energy]= textread( 'TotalEnergy.txt',...
    repmat('%f ',[1,2]), 'headerlines',1);

plot(time, energy, 'lineWidth', 3);