function [ output_args ] = AddParticleStream(num, x0, y0, PartAng, Type, Ep, Seper, shape, velocityChange, velocityChangeFactor)
global AtomSpacing x y AtomType Vx Vy Mass0 Mass1 nAtoms

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end
%Shape == 0 - linear particle stream
%Shape ==1 - rect particle stream
%Shape ==2 - V
if shape == 0
    for p = 0:num - 1
        nAtoms = nAtoms + 1;
        x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng);
        y(nAtoms) = y0 * AtomSpacing - Seper * p * AtomSpacing * sin(PartAng);
        AtomType(nAtoms) = Type;
    end
elseif shape ==1
    for p=0:num -1
        if mod(p,2)==0
            nAtoms = nAtoms + 1;
            x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng) + AtomSpacing;
            y(nAtoms) = y0 * AtomSpacing - Seper * p * AtomSpacing * sin(PartAng);
            AtomType(nAtoms) = Type;   
        else
            nAtoms = nAtoms + 1;
            x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng) - AtomSpacing;
            y(nAtoms) = y0 * AtomSpacing - Seper * (p-1) * AtomSpacing * sin(PartAng);
            AtomType(nAtoms) = Type;   
        end
    end
elseif shape==2
    for p=0:num -1
        if(p==0)
            nAtoms = nAtoms + 1;
            x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng);
            y(nAtoms) = y0 * AtomSpacing - Seper * p * AtomSpacing * sin(PartAng);
            AtomType(nAtoms) = Type;
        else
            if mod(p,2)~=0
                nAtoms = nAtoms + 1;
                x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng) + p*AtomSpacing;
                y(nAtoms) = y0 * AtomSpacing - Seper * p * AtomSpacing * sin(PartAng);
                AtomType(nAtoms) = Type;   
            else
                nAtoms = nAtoms + 1;
                x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng) - (p-1)*AtomSpacing;
                y(nAtoms) = y0 * AtomSpacing - Seper * (p-1) * AtomSpacing * sin(PartAng);
                AtomType(nAtoms) = Type;   
            end
        end        
    end
end
V = sqrt(2 * Ep / Mass);

for p = 1:num
    if velocityChange ==0 %Constant
        Vx(nAtoms - num + p) = V * cos(PartAng)/velocityChangeFactor;
        Vy(nAtoms - num + p) = V * sin(PartAng)/velocityChangeFactor;
    elseif velocityChange ==1 %Linear increase
        Vx(nAtoms - num + p) = V * cos(PartAng)*p/velocityChangeFactor;
        Vy(nAtoms - num + p) = V * sin(PartAng)*p/velocityChangeFactor;
    elseif velocityChange ==2 %Linear Decrease
        Vx(nAtoms - num + p) = V * cos(PartAng)*(1/p)/velocityChangeFactor;
        Vy(nAtoms - num + p) = V * sin(PartAng)*(1/p)/velocityChangeFactor;  
    elseif velocityChange ==3 %Quadratic Increase
        Vx(nAtoms - num + p) = V * cos(PartAng)*(p^2)/velocityChangeFactor;
        Vy(nAtoms - num + p) = V * sin(PartAng)*(p^2)/velocityChangeFactor;
    elseif velocityChange ==4 %Quadratoc Decrease
        Vx(nAtoms - num + p) = V * cos(PartAng)*(1/(p^2))/velocityChangeFactor;
        Vy(nAtoms - num + p) = V * sin(PartAng)*(1/(p^2))/velocityChangeFactor;
    elseif velocityChange ==5 %Sinusoidal
        Vx(nAtoms - num + p) = V * cos(PartAng)*abs(sin(p))/velocityChangeFactor;
        Vy(nAtoms - num + p) = V * sin(PartAng)*abs(sin(p))/velocityChangeFactor;
    end
end

end
