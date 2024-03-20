use std::fs;
use std::path::Path;
///
/// 测试方程文件的可靠性
///
#[test]
fn test_json() {
    let res_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("res");
    let res_dir = fs::read_dir(res_path);
    let mut vec_file: Vec<String> = vec![];
    for entry in res_dir.expect("no res") {
        if let Ok(entry) = entry {
            if entry.file_name().to_string_lossy().contains(".json") {
                vec_file.push(entry.file_name().to_string_lossy().to_string());
            }
        }
    }
    println!("{:?}", vec_file);
}
