% script:  implement_bessel_digital_filter.m
% descrip: Used to prove-in the C++ implementation of digital bessel filters

% Defs 
Norders = 3;

% filters
bessel_poly{1} = [1 1];
bessel_poly{2} = [1 3 3];
bessel_poly{3} = [1 6 15 15];
bessel_poly{4} = [1 10 45 105 105];
bessel_poly{5} = [1 15 105 420 945 945];

xi_3db = [1 1.361 1.756 2.114 2.427];
xi_scale = [1 1 1 1 1];

% axes
T = 1;
Neff = 100;

% Neff = 1, xi_scale = 1 pole and zero locations
for order = 3:Norders,
    
    % analog poles and residues
    [res{order}, p{order}, k] = residue(max(bessel_poly{order}), bessel_poly{order});

    % digital poles
    zpole{order} = exp(p{order} * T ./ (xi_3db(order) * xi_scale(order) * (Neff / xi_3db(order))));
    
    switch order
        case 1,
            ycoef = [1 -zpole{order}];
        case 2,
            ycoef = conv([1 -zpole{order}(1)], [1 -zpole{order}(2)]);
        case 3,
            ycoef = conv(conv([1 -zpole{order}(1)], [1 -zpole{order}(2)]), [1 -zpole{order}(3)]);
        case 4,
            ycoef = conv(conv(conv([1 -zpole{order}(1)], [1 -zpole{order}(2)]), [1 -zpole{order}(3)]), [1 -zpole{order}(4)]);
        case 5,
            ycoef = conv(conv(conv(conv([1 -zpole{order}(1)], [1 -zpole{order}(2)]), [1 -zpole{order}(3)]), [1 -zpole{order}(4)]), [1 -zpole{order}(5)]);
    end
        
    ycoef = -real(ycoef(2:end));
    s = ['order ' num2str(order) ' Ycoef: '];
    for j = 1: order,
        s = [s num2str(ycoef(j), 13) '  '];
    end
    disp(s)
    
    % zeros
    switch order
        case 1,
            xcoef = res{order} * [1 +zpole{order}];
        case 2,
            xcoef = ...
                res{order}(1) * conv([1                 ], [1 -zpole{order}(2)]) + ...
                res{order}(2) * conv([1 -zpole{order}(1)], [1                 ]);
        case 3,
            xcoef = ...
                res{order}(1) * conv([1                 ], conv([1 -zpole{order}(2)], [1 -zpole{order}(3)])) + ...
                res{order}(2) * conv([1 -zpole{order}(1)], conv([1                 ], [1 -zpole{order}(3)])) + ...
                res{order}(3) * conv([1 -zpole{order}(1)], conv([1 -zpole{order}(2)], [1                 ]));
        case 4,
            xcoef = ...
                res{order}(1) * conv([1                 ], conv([1 -zpole{order}(2)], conv([1 -zpole{order}(3)], [1 -zpole{order}(4)]))) + ...
                res{order}(2) * conv([1 -zpole{order}(1)], conv([1                 ], conv([1 -zpole{order}(3)], [1 -zpole{order}(4)]))) + ...
                res{order}(3) * conv([1 -zpole{order}(1)], conv([1 -zpole{order}(2)], conv([1                 ], [1 -zpole{order}(4)]))) + ...
                res{order}(4) * conv([1 -zpole{order}(1)], conv([1 -zpole{order}(2)], conv([1 -zpole{order}(3)], [1                 ])));
        case 5,
            xcoef = ...
                res{order}(1) * conv([1                 ], conv([1 -zpole{order}(2)], conv([1 -zpole{order}(3)], conv([1 -zpole{order}(4)], [1 -zpole{order}(5)])))) + ...
                res{order}(2) * conv([1 -zpole{order}(1)], conv([1                 ], conv([1 -zpole{order}(3)], conv([1 -zpole{order}(4)], [1 -zpole{order}(5)])))) + ...
                res{order}(3) * conv([1 -zpole{order}(1)], conv([1 -zpole{order}(2)], conv([1                 ], conv([1 -zpole{order}(4)], [1 -zpole{order}(5)])))) + ...
                res{order}(4) * conv([1 -zpole{order}(1)], conv([1 -zpole{order}(2)], conv([1 -zpole{order}(3)], conv([1                 ], [1 -zpole{order}(5)])))) + ...
                res{order}(5) * conv([1 -zpole{order}(1)], conv([1 -zpole{order}(2)], conv([1 -zpole{order}(3)], conv([1 -zpole{order}(4)], [1                 ]))));
    end
    
    gain = T / (xi_3db(order) * xi_scale(order) * Neff) * real(res{order}.' * (1./(1 - zpole{order})));
    
    xcoef = T / (xi_3db(order) * xi_scale(order) * Neff) * real(xcoef);
    s = ['order ' num2str(order) ' Xcoef: '];
    for j = 1: order,
        s = [s num2str(xcoef(j), 13) '  '];
    end
    disp(s)
    disp(' ')
    
end






