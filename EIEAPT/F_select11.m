% The selection function in RVEA
function [Selection] = F_select11(FunctionValue, V, theta0)

[N M] = size(FunctionValue);
VN = size(V, 1);

Zmin = min(FunctionValue,[],1);

%Translation
FunctionValue = (FunctionValue - repmat(Zmin, [size(FunctionValue,1) 1]));

%Solutions associattion to reference vectors
clear class;
uFunctionValue = FunctionValue./repmat(sqrt(sum(FunctionValue.^2,2)), [1 M]);
cosine = uFunctionValue*V'; %calculate the cosine values between each solution and each vector
acosine = acos(cosine);
[maxc maxcidx] = max(cosine, [], 2);
class = struct('c', cell(1,VN)); %classification
for i = 1:N
    class(maxcidx(i)).c = [class(maxcidx(i)).c, i];
end;



Selection = [];

for k = 1:VN
    if(~isempty(class(k).c))
        sub = class(k).c;
        subFunctionValue = FunctionValue(sub,:);
        
        
        D=calculate_pbi_matrix_optimized(subFunctionValue,V(k,:),theta0(k));
        %PBI calculation
%         subacosine = acosine(sub, k);
%         subacosine = subacosine/refV(k); % angle normalization
%         D1 = sqrt(sum(subFunctionValue.^2,2)); % Euclidean distance from solution to the ideal point
%         D2=D1*sqrt(1-cosine(:,k).^2);
%         D = D1+theta(k)*D2; % APD
%          D = D1.*(1 + (theta0)*(subacosine)); % APD
        
        [mind mindidx] = min(D);
        Selection = [Selection; sub(mindidx)];
    end;
end;

end
function pbi_values = calculate_pbi_matrix_optimized(solutions, reference_direction, theta)
    % solutions: �����ÿһ�д���һ�����Ŀ��ֵ����
    % reference_direction: �ο���������
    % theta: �ͷ�����

    % ��ȡ���������Ŀ�������
    [num_solutions, num_objectives] = size(solutions);
    
    % �������������������Ҫ��ʽ���壬��Ϊ��ȥ������û��Ч��
    % ideal_point = zeros(1, num_objectives);
    
    % ����ο����������ķ���
    norm_reference = norm(reference_direction);
    
    % ���������������
    d = solutions; % solutions - ideal_point����ȥ������û��Ч��
    
    % ������������ڲο������ϵ�ͶӰ���Ⱦ���
    projection_length = (d * reference_direction') / norm_reference;
    
    % ����ͶӰ��������
    projection_vector = (projection_length / norm_reference) * reference_direction;
    
    % ���㴹ֱ�ڲο��������������
    orthogonal_vector = d - projection_vector;
    
    % ����ÿ�����PBIֵ
    pbi_values = sqrt(sum(projection_vector.^2, 2)) + theta * sqrt(sum(orthogonal_vector.^2, 2));
end
