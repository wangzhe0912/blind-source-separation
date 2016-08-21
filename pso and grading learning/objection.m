function [ f ] = objection ( newDi,Di )
%计算目标函数
f=(newDi-Di)/Di;
f=max(f);
end

