vector[] ode_explicit_trapezoidal(vector initial_state,
                                    real initial_time,
                                    real[] times,
                                    real beta,
                                    real omega,
                                    real dI,
                                    real dT) {
    real h;
    vector[size(initial_state)] dstate_dt_initial_time;
    vector[size(initial_state)] dstate_dt_tidx;
    vector[size(initial_state)] k;
    vector[size(initial_state)] state_estimate[size(times)];

    h = times[1] - initial_time;
    dstate_dt_initial_time = sird(initial_time, initial_state, beta, omega, dI, dT);
    k = h * dstate_dt_initial_time;
    state_estimate[1,] = initial_state + h * (dstate_dt_initial_time + sird(times[1], initial_state + k, beta, omega, dI, dT)) / 2;

    for (tidx in 1:size(times)-1) {
      h = (times[tidx+1] - times[tidx]);
      dstate_dt_tidx = sird(times[tidx], state_estimate[tidx], beta, omega, dI, dT);
      k = h * dstate_dt_tidx;
      state_estimate[tidx+1,] = state_estimate[tidx,] + h * (dstate_dt_tidx + sird(times[tidx+1], state_estimate[tidx,] + k, beta, omega, dI, dT))/2;
    }

    return state_estimate;
  }
