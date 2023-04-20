% A collection of utility functions




  function check_args(arg_vector)
  % check_args
  %   checks if the vector of variable names exists.  Used to validate args exist
  %   reports missing args and raises an error if any args are missing
  % ARGS:
  %   arg_vector: vector, variables to check for existence.
    E_missing_arg = false;
    for arg = arg_vector 
      if ~exist(arg,'var')
        fprintf('MISSING ARG: %s\n', arg)
        E_missing_arg = true;
      end
    end
    if E_missing_arg
      error('ERROR: Missing Arguments!')
    end
  end
