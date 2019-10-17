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
% Model name = pancreas_min_1
%
%


function main()
%Initial conditions vector
	x0=zeros(6,1);
	x0(1) = 5.0;
	x0(2) = 0.8;
	x0(3) = 6.0E-5;
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
	compartment_Vext=0.1;
% Parameter:   id =  Cext_glc, name = glucose concentration (Vext)
% Parameter:   id =  Cext_lac, name = lactate concentration (Vext)
% Parameter:   id =  Cext_ins, name = insulin concentration (Vext)
% Parameter:   id =  Cext_cpep, name = c-peptide concentration (Vext)
% Parameter:   id =  Cpa_glc, name = glucose concentration (Vpa)
% Parameter:   id =  Cpa_lac, name = lactate concentration (Vpa)
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
% assignmentRule: variable = Cext_glc
	global_par_Cext_glc=x(1)/compartment_Vext;
% assignmentRule: variable = Cext_lac
	global_par_Cext_lac=x(2)/compartment_Vext;
% assignmentRule: variable = Cext_ins
	global_par_Cext_ins=x(3)/compartment_Vext;
% assignmentRule: variable = Cext_cpep
	global_par_Cext_cpep=x(4)/compartment_Vext;
% assignmentRule: variable = Cpa_glc
	global_par_Cpa_glc=x(5)/compartment_Vpa;
% assignmentRule: variable = Cpa_lac
	global_par_Cpa_lac=x(6)/compartment_Vpa;

% Reaction: id = GLCIM, name = glucose import
	reaction_GLCIM=compartment_Vpa*global_par_GLCIM_Vmax/global_par_GLCIM_Km*(global_par_Cext_glc-global_par_Cpa_glc)/(1+global_par_Cext_glc/global_par_GLCIM_Km+global_par_Cpa_glc/global_par_GLCIM_Km);

% Reaction: id = LACEX, name = lactate export
	reaction_LACEX=compartment_Vpa*global_par_LACEX_Vmax/global_par_LACEX_Km*(global_par_Cpa_lac-global_par_Cext_lac)/(1+global_par_Cext_lac/global_par_LACEX_Km+global_par_Cpa_lac/global_par_LACEX_Km);

% Reaction: id = GLC2LAC, name = glycolysis
	reaction_GLC2LAC=compartment_Vpa*global_par_GLC2LAC_Vmax*global_par_Cpa_glc/(global_par_Cpa_glc+global_par_GLC2LAC_Km);

% Reaction: id = IRS, name = IRS insulin secretion
	reaction_IRS=compartment_Vpa*global_par_IRS_Vmax*global_par_Cpa_glc^global_par_IRS_n_glc/(global_par_Cpa_glc^global_par_IRS_n_glc+global_par_IRS_Km_glc^global_par_IRS_n_glc);

	xdot=zeros(6,1);
	
% Species:   id = Aext_glc, name = glucose
%WARNING speciesID: Aext_glc, constant= false  , boundaryCondition = true but is not involved in assignmentRule, rateRule or events !
	xdot(1) = 0.0;
	
% Species:   id = Aext_lac, name = lactate, affected by kineticLaw
	xdot(2) = ( 1.0 * reaction_LACEX);
	
% Species:   id = Aext_ins, name = insulin, affected by kineticLaw
	xdot(3) = ( 1.0 * reaction_IRS);
	
% Species:   id = Aext_cpep, name = c-peptide, affected by kineticLaw
	xdot(4) = ( 1.0 * reaction_IRS);
	
% Species:   id = Apa_glc, name = glucose, affected by kineticLaw
	xdot(5) = ( 1.0 * reaction_GLCIM) + (-1.0 * reaction_GLC2LAC);
	
% Species:   id = Apa_lac, name = lactate, affected by kineticLaw
	xdot(6) = (-1.0 * reaction_LACEX) + ( 2.0 * reaction_GLC2LAC);
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


