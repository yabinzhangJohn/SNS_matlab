
% example of how to call gbvs with default params

% img = imread('samplepics/4.jpg');
% img = All_original_im{1};
% img = imread('D:\ImageRetargeting\MIT dataset\ArtRoom\ArtRoom.png');


PATH_NAME = {'ArtRoom'; 'BedRoom'; 'Brasserie_L_Aficion'; ...
    'DKNYgirl'; 'Deck'; 'Fatem'; 'Johanneskirche'; 'Lotus'; 'Marblehead_Mass';...
    'Perissa_Santorini'; 'Sanfrancisco'; 'SetAngle'; 'Umdan'; 'Unazukin';...
    'Woman'; 'bicycle2'; 'brick_house'; 'buddha'; 'butterfly'; 'car1'; 'car';...
    'child'; 'face'; 'family'; 'foliage'; 'getty'; 'girls'; 'glasses'; 'greek_wine';...
    'jon'; 'mnm'; 'obama'; 'painting2'; 'surfers'; 'tajmahal'; 'tower'; 'volleyball'};
PATH_ROOT = 'D:\ImageRetargeting\MIT dataset\';
for set_num = 1:37
    path = [PATH_ROOT PATH_NAME{set_num} '\'];
    file = dir([path,'*.png']);
    imgList{set_num} = imread([path file(1).name]);
end

% for set_num = 1:35
%     prefix = 'D:\ImageRetargeting\LinCW\data\';
%     type = {'mo'; 'scaling'; 'seam'; 'shift'; 'warp'; 'src'};
%     original_path = [prefix type{6} '\'];
%     % Extension name of the input file
%     ext_type = '*.jpg';
%     img_list = dir([original_path ext_type]);   %src
%     original_im = imread([original_path img_list(set_num).name]);
%     imgList{set_num} = original_im;
% end

for i = 4
    img = imgList{i};
    out_gbvs = gbvs(img);
    out_itti = ittikochmap(img);

    figure;
    subplot(2,3,1);
    imshow(img);
    title('Original Image');

    subplot(2,3,2);
    show_imgnmap(img,out_gbvs);
    title('GBVS map overlayed');

    subplot(2,3,3);
    show_imgnmap(img,out_itti);
    title('Itti/Koch map overlayed');


    subplot(2,3,5);
    imshow( out_gbvs.master_map_resized );
    title('GBVS map');

    subplot(2,3,6);
    imshow(out_itti.master_map_resized);
    title('Itti/Koch map');
end
