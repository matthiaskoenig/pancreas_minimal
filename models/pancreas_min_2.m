% This file works with MATLAB and is automatically generated with 
% the System Biology Format Converter (http://sbfc.sourceforge.net/)
% from an SBML file. 
% To run this file with Octave you must edit the comments providing
% the definition of the ode solver and the signature for the 
% xdot function.
%
% The conversion system has the following limitations:
%  - You may have to re order some reactions and Assignment Rules definition
%  - Delays are not taken into account
%  - You should change the lsode parameters (start, end, steps) to get better results
%

%
% Model name = pancreas_min_2
%
%


function main()
%Initial conditions vector
	x0=zeros(6,1);
	x0(1) = 5.0;
	x0(2) = 0.8;
	x0(3) = 6.0E-8;
	x0(4) = 0.0;
	x0(5) = 5.0;
	x0(6) = 0.8;


% Depending on whether you are using Octave or Matlab,
% you should comment / uncomment one of the following blocks.
% This should also be done for the definition of the function f below.
% Start Matlab code
	tspan=[0:0.01:100];
	opts = odeset('AbsTol',1e-3);
	[t,x]=ode23tb(@f,tspan,x0,opts);
% End Matlab code

% Start Octave code
%	t=linspace(0,100,100);
%	x=lsode('f',x0,t);
% End Octave code


	plot(t,x);
end



% Depending on whether you are using Octave or Matlab,
% you should comment / uncomment one of the following blocks.
% This should also be done for the definition of the function f below.
% Start Matlab code
function xdot=f(t,x)
% End Matlab code

% Start Octave code
%function xdot=f(x,t)
% End Octave code

% Compartment: id = Vpa, name = pancreas tissue, constant
	compartment_Vpa=0.5;
% Compartment: id = Vext, name = pancreas blood, constant
	compartment_Vext=5.0;
% Parameter:   id =  GLCIM_Vmax, name = Glucose import
	global_par_GLCIM_Vmax=100.0;
% Parameter:   id =  GLCIM_Km, name = GLCIM_Km
	global_par_GLCIM_Km=1.0;
% Parameter:   id =  LACEX_Vmax, name = Lactate import
	global_par_LACEX_Vmax=100.0;
% Parameter:   id =  LACEX_Km, name = LACEX_Km
	global_par_LACEX_Km=0.5;
% Parameter:   id =  GLC2LAC_Vmax, name = Glucose utilization
	global_par_GLC2LAC_Vmax=0.1;
% Parameter:   id =  GLC2LAC_Km, name = GLC2LAC_Km
	global_par_GLC2LAC_Km=4.5;
% Parameter:   id =  IRS_Vmax, name = Insulin secretion
	global_par_IRS_Vmax=1.6E-6;
% Parameter:   id =  IRS_n_glc, name = IRS_n_glc
	global_par_IRS_n_glc=4.0;
% Parameter:   id =  IRS_Km_glc, name = IRS_Km_glc
	global_par_IRS_Km_glc=7.0;

% Reaction: id = GLCIM, name = glucose import
	reaction_GLCIM=compartment_Vpa*global_par_GLCIM_Vmax/global_par_GLCIM_Km*(x(1)-x(5))/(1+x(1)/global_par_GLCIM_Km+x(5)/global_par_GLCIM_Km);

% Reaction: id = LACEX, name = lactate export
	reaction_LACEX=compartment_Vpa*global_par_LACEX_Vmax/global_par_LACEX_Km*(x(6)-x(2))/(1+x(2)/global_par_LACEX_Km+x(6)/global_par_LACEX_Km);

% Reaction: id = GLC2LAC, name = glycolysis
	reaction_GLC2LAC=compartment_Vpa*global_par_GLC2LAC_Vmax*x(5)/(x(5)+global_par_GLC2LAC_Km);

% Reaction: id = IRS, name = IRS insulin secretion
	reaction_IRS=compartment_Vpa*global_par_IRS_Vmax*x(5)^global_par_IRS_n_glc/(x(5)^global_par_IRS_n_glc+global_par_IRS_Km_glc^global_par_IRS_n_glc);

	xdot=zeros(6,1);
	
% Species:   id = Cext_glc, name = glucose
%WARNING speciesID: Cext_glc, constant= false  , boundaryCondition = true but is not involved in assignmentRule, rateRule or events !
	xdot(1) = 0.0;
	
% Species:   id = Cext_lac, name = lactate, affected by kineticLaw
	xdot(2) = (1/(compartment_Vext))*(( 1.0 * reaction_LACEX));
	
% Species:   id = Cext_ins, name = insulin, affected by kineticLaw
	xdot(3) = (1/(compartment_Vext))*(( 1.0 * reaction_IRS));
	
% Species:   id = Cext_cpep, name = c-peptide, affected by kineticLaw
	xdot(4) = (1/(compartment_Vext))*(( 1.0 * reaction_IRS));
	
% Species:   id = Cpa_glc, name = glucose, affected by kineticLaw
	xdot(5) = (1/(compartment_Vpa))*(( 1.0 * reaction_GLCIM) + (-1.0 * reaction_GLC2LAC));
	
% Species:   id = Cpa_lac, name = lactate, affected by kineticLaw
	xdot(6) = (1/(compartment_Vpa))*((-1.0 * reaction_LACEX) + ( 2.0 * reaction_GLC2LAC));
end

% adding few functions representing operators used in SBML but not present directly 
% in either matlab or octave. 
function z=pow(x,y),z=x^y;end
function z=root(x,y),z=y^(1/x);end
function z = piecewise(varargin)
	numArgs = nargin;
	result = 0;
	foundResult = 0;
	for k=1:2: numArgs-1
		if varargin{k+1} == 1
			result = varargin{k};
			foundResult = 1;
			break;
		end
	end
	if foundResult == 0
		result = varargin{numArgs};
	end
	z = result;
end


