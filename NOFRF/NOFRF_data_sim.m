function Y = NOFRF_data_sim(ode_func,tspan,n_A,A,Y,y0)

Ts = mean(abs(diff(tspan)));
for i = 1:n_A
    clear y1
    Amp_1 = A(i);
    y1 = ode4(@(t,y) ode_func(t,y,(Amp_1)),tspan,y0);
    % y1 = RK4_solver(@(t,y) ode_func(t,y,(Amp_1)) , [tspan(1), tspan(end)], y0, Ts);
    Y(:,i) = y1(:,1);
end

end


