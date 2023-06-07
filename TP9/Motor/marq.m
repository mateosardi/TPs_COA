function [W1,W2,PI_vector,iteration,lambda]=marq(NetDef,W1,W2,PHI,Y,trparms)
%  MARQ
%  ----
%              Train a two layer neural network with the Levenberg-Marquardt
%              method.
%
%          If desired, it is possible to use regularization by
%          weight decay. Also pruned (ie. not fully connected) networks can
%          be trained.
%
%          Given a set of corresponding input-output pairs and an initial
%          network,
%          [W1,W2,critvec,iteration,lambda]=marq(NetDef,W1,W2,PHI,Y,trparms)
%          trains the network with the Levenberg-Marquardt method.
% 
%          The activation functions can be either linear or tanh. The
%          network architecture is defined by the matrix 'NetDef' which
%          has two rows. The first row specifies the hidden layer and the
%          second row specifies the output layer.
% 
%          E.g.:    NetDef = ['LHHHH' 
%                             'LL---']
%
%          (L = Linear, H = tanh)
%
%          A weight is pruned by setting it to zero.
%
%          The Marquardt method is described in:
%          K. Madsen: 'Optimering' (Haefte 38), IMM, DTU, 1991
%  
%          Notice that the bias is included as the last column in the weight
%          matrices.
%
% 
%  INPUT:
%  NetDef: Network definition 
%  W1    : Input-to-hidden layer weights. The matrix dimension is
%          dim(W1) = [(# of hidden units) * (inputs + 1)] (the 1 is due to the bias)
%  W2    : hidden-to-output layer weights.
%           dim(W2) = [(outputs)  *  (# of hidden units + 1)]
%  PHI   : Input vector. dim(PHI) = [(inputs)  *  (# of data)]
%  Y     : Output data. dim(Y) = [(outputs)  * (# of data)]
%  trparms : Vector containing parameters associated with the training
%            trparms = [max_iter stop_crit lambda D]
%             max_iter  : max # of iterations.
%             stop_crit : Stop training if criterion is below this value
%             lambda    : Initial Levenberg-Marquardt parameter
%             D         : Row vector containing the weight decay parameters.
%                            If D has one element a scalar weight decay will be
%                            used. If D has two elements the first element will
%                            be used as weight decay for the hidden-to-output
%                            layer while second will be used for the input-to
%                            hidden layer weights. For individual weight decays,
%                            D must contain as many elements as there are
%                            weights in the network.
%
%           Default values are (obtained if left out): trparms = [500 0 1 0] 
% 
% 
%  OUTPUT:
%  W1, W2   : Weight matrices after training
%  critvec:   Vector containing the criterion evaluated at each iteration
%  iteration: # of iterations
%  lambda   : The final value of lambda. Relevant only if retraining is desired
% 
%  Programmed by : Magnus Norgaard, IAU/IMM
%  LastEditDate  : July 16, 1996


%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
[outputs,N] = size(Y);                  % # of outputs and # of data
[hidden,inputs] = size(W1);             % # of hidden units 
inputs=inputs-1;                        % # of inputs
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = [zeros(hidden,N);ones(1,N)]; % Hidden layer outputs
y2       = zeros(outputs,N);            % Network output
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A useful vector!
index2 = (0:N-1)*outputs;               % Yet another useful vector
iteration = 1;                          % Counter variable
dw       = 1;                           % Flag telling that the weights are new
PHI      = [PHI;ones(1,N)];             % Augment PHI with a row containg ones
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters = parameters1 + parameters2; % Total # of weights
PSI      = zeros(parameters,outputs*N); % Deriv. of each output w.r.t. each weight
ones_h   = ones(hidden+1,1);            % A vector of ones
ones_i   = ones(inputs+1,1);            % Another vector of ones
                                        % Parameter vector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3   = 1:(reduced+1):(reduced^2);   % A third useful vector
if ~exist('trparms')                    % Default training parameters
  max_iter  = 500;
  stop_crit = 0;
  lambda    = 1;
  D         = 0;
else                                    % User specified values
  max_iter  = trparms(1);
  stop_crit = trparms(2);
  lambda    = trparms(3);
  if length(trparms)==4,                % Scalar weight decay parameter
    D = trparms(4*ones(1,reduced))';      
  elseif length(trparms)==5,            % Two weight decay parameters
    D = trparms([4*ones(1,parameters2) 5*ones(1,parameters1)])';
    D = D(theta_index);
  elseif length(trparms)>5,             % Individual weight decay
    D = trparms(4:length(trparms))';
  end
end
PI_vector = zeros(max_iter,1);          % A vector containing the accumulated SSE


%----------------------------------------------------------------------------------
%--------------                   TRAIN NETWORK                       -------------
%----------------------------------------------------------------------------------
clc;
c=fix(clock);
fprintf('Network training started at %2i.%2i.%2i\n\n',c(4),c(5),c(6));


% >>>>>>>>>>>>>>>>>>>>>  COMPUTE NETWORK OUTPUT  y2(theta)   <<<<<<<<<<<<<<<<<<<<<<
h1 = W1*PHI;  
y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
y1(L_hidden,:) = h1(L_hidden,:);    

h2 = W2*y1;
y2(H_output,:) = pmntanh(h2(H_output,:));
y2(L_output,:) = h2(L_output,:);

E        = Y - y2;                      % Training error
E_vector = E(:);                        % Reshape E into a long vector
SSE      = E_vector'*E_vector;          % Sum of squared errors (SSE)
PI       = (SSE+theta_red'*(D.*theta_red))/(2*N); % Performance index

while iteration<=max_iter    
if dw==1,
% >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
% (The derivative of each network output (y2) with respect to each weight)

    % ==========   Elements corresponding to the linear output units   ============
    for i = L_output'
      index1 = (i-1) * (hidden + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      PSI(index1:index1+hidden,index2+i) = y1;
      % ---------------------------------------------------------------------
 
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        PSI(index(j):index(j)+inputs,index2+i) = W2(i,j)*PHI;
      end

      for j = H_hidden',
        tmp = W2(i,j)*(1-y1(j,:).*y1(j,:)); 
        PSI(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:).*PHI;
      end 
      % ---------------------------------------------------------------------    
    end

    % ============  Elements corresponding to the tanh output units   =============
    for i = H_output',
      index1 = (i-1) * (hidden + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      tmp = 1 - y2(i,:).*y2(i,:);
      PSI(index1:index1+hidden,index2+i) = y1.*tmp(ones_h,:);
      % ---------------------------------------------------------------------
     
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        tmp = W2(i,j)*(1-y2(i,:).*y2(i,:));
        PSI(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:).*PHI;
      end
      
      for j = H_hidden',
        tmp  = W2(i,j)*(1-y1(j,:).*y1(j,:));
        tmp2 = (1-y2(i,:).*y2(i,:));
        PSI(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:)...
                                                  .*tmp2(ones_i,:).*PHI;
      end
      % ---------------------------------------------------------------------
    end
    PSI_red = PSI(theta_index,:);
    
    % -- Gradient --
    G = PSI_red*E_vector-D.*theta_red;

    % -- Means square error part Hessian  --
    R = PSI_red*PSI_red';
    
    dw = 0;
  end
  
   
% >>>>>>>>>>>>>>>>>>>>>>>>>>>        COMPUTE h_k        <<<<<<<<<<<<<<<<<<<<<<<<<<<
  % -- Hessian  --
  H = R;
  H(index3) = H(index3)'+lambda+D;                  % Add diagonal matrix

  % -- Search direction --
  h = H\G;                                          % Solve for search direction

  % -- Compute 'apriori' iterate --
  theta_red_new = theta_red + h;                    % Update parameter vector
  theta(theta_index) = theta_red_new;

  % -- Put the parameters back into the weight matrices --
  W1_new = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
  W2_new = reshape(theta(1:parameters2),hidden+1,outputs)';

    
% >>>>>>>>>>>>>>>>>>>>   COMPUTE NETWORK OUTPUT  y2(theta+h)   <<<<<<<<<<<<<<<<<<<<
  h1 = W1_new*PHI;  
  y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
  y1(L_hidden,:) = h1(L_hidden,:);
    
  h2 = W2_new*y1;
  y2(H_output,:) = pmntanh(h2(H_output,:));
  y2(L_output,:) = h2(L_output,:);

  E_new        = Y - y2;                 % Training error
  E_new_vector = E_new(:);               % Reshape E into a long vector
  SSE_new  = E_new_vector'*E_new_vector; % Sum of squared errors (SSE)
  PI_new   = (SSE_new + theta_red_new'*(D.*theta_red_new))/(2*N); % PI


% >>>>>>>>>>>>>>>>>>>>>>>>>>>       UPDATE  lambda     <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  L = h'*G + h'*(h.*(D+lambda));

  % Decrease lambda if SSE has fallen 'sufficiently'
  if 2*N*(PI - PI_new) > (0.75*L),
    lambda = lambda/2;
  
  % Increase lambda if SSE has grown 'sufficiently'
  elseif 2*N*(PI-PI_new) <= (0.25*L),
    lambda = 2*lambda;
  end


% >>>>>>>>>>>>>>>>>>>>       UPDATES FOR NEXT ITERATION        <<<<<<<<<<<<<<<<<<<<
  % Update only if criterion has decreased
  if PI_new < PI,                      
    W1 = W1_new;
    W2 = W2_new;
    theta_red = theta_red_new;
    E_vector = E_new_vector;
    PI = PI_new;
    dw = 1;
    iteration = iteration + 1;
    PI_vector(iteration-1) = PI;                             % Collect PI in vector
    fprintf('iteration # %i   PI = %4.3e\r',iteration-1,PI); % Print on-line inform
  end

  % Check if stop condition is fulfilled
  if (PI < stop_crit) | (lambda>1e7), break, end             
end
%----------------------------------------------------------------------------------
%--------------              END OF NETWORK TRAINING                  -------------
%----------------------------------------------------------------------------------
PI_vector = PI_vector(1:iteration-1);
c=fix(clock);
fprintf('\n\nNetwork training ended at %2i.%2i.%2i\n',c(4),c(5),c(6));
