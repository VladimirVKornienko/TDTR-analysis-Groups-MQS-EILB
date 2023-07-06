function [t_interior,f_interior] = extract_interior_V4(t_data,f_data,t_min,t_max)
%EXTRACTS all points in range tmin<= t_data <= tmax, and corresponding
%points of f (of same length)
    Ndata=length(t_data);
    Nmin=Ndata;
    for i=Ndata:-1:1
        if t_data(i)>=t_min
            Nmin=i;
        end
    end
    Nmax=Nmin;
    for i=Nmin:Ndata
        if t_data(i)<=t_max
            Nmax=i;
        end
    end
    t_interior=t_data(Nmin:Nmax,1);
    f_interior=f_data(Nmin:Nmax,1);
end

