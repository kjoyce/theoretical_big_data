classdef FinSupFun % Finite Support Function
    properties (SetAccess = private)
        l =  Inf;  
        r = -Inf;
        f = [];   % Function w support [l,r]
    end
    
    methods
        function c = FinSupFun(F,L) % Constructor
            if nargin == 0  % Empty function
                c.l =  Inf;  
                c.r = -Inf;
                c.f =  [];   % Function w support [l,r]
            elseif nargin == 1   % Symmetric (even) function
                c.r = size(F,2)-1; 
                c.l = -c.r; 
                c.f = [fliplr(F) F(2:c.r+1)];
%                 c.l = 0; 
%                 c.r = size(F,2)-1; 
%                 c.f = F;
            elseif nargin == 2
                c.l = L; 
                c.r = L+size(F,2)-1; 
                c.f = F;
            end
        end
        
        function c = plus(a,b)  % Addition
            L = min(a.l, b.l);
            R = max(a.r, b.r);
            F = zeros(1, R-L+1);
            F((a.l:a.r)-L+1) = a.f;
            F((b.l:b.r)-L+1) = F((b.l:b.r)-L+1) + b.f;
            c = FinSupFun(F, L);
        end
        
        function c = mtimes(a,b)  % * Convolution Full
            c = FinSupFun(conv(a.f,b.f), a.l + b.l);
        end
        
        function c = times(a,b)  % .* Convolution Internal
            c = FinSupFun(conv(b.f,a.f,'valid'), b.l + a.r);
        end
        
        function b = ctranspose(a)  % ' Reverse (transposed)
            b = FinSupFun(fliplr(a.f), -a.r);
        end

	function c = restricted_to(a,l,r) % Restrict support to [l,r], if [l,r] is bigger than [a.l,a.r], then pad with zeroes.
	    L = max(l,a.l);	% Left endpoint of restricted interval 
	    R = min(r,a.r);	% Right endpoint of restricted interval 
	    start_idx = L - a.l + 1; % How many points to cut from the left
	    len = R - L; % How far from the left to include
	    b = FinSupFun(a.f(start_idx:(start_idx+len)),L); % Define a on restricted interval
	    c = b + FinSupFun(zeros(1,r-l+1),l); % Now, pad with zeros using logic in +
	end

	function c = mldivide(a,b) % \ De-convolution by constructing toeplitz matrix. a must be symmetric 
	  n = length(b.f);
	  toeplitz_row = [a.f((a.r+1):end), zeros(1, n-(a.r))];% This needs to be from the center of p and padded with zeros
	  if length(toeplitz_row) > length(b.f)
	    toeplitz_row = toeplitz_row(1:length(b.f));
	  end
	  A = toeplitz(toeplitz_row);
	  c = FinSupFun((A\b.f')',b.l); 
	end
    end
end % classdef

